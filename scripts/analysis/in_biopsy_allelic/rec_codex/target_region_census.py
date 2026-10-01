"""Outcome-blind Currin target geometry census on development folds 1--4."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
LABELS = ROOT / "GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz"
PEAKS = ROOT / "GWAS/finemapping/data/seqfunc_external/currin2025_caqtl_v1/supplementalData1_liver_ATAC_peaks.bed.gz"
FAI = Path("/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa.fai")
FIELDS = ["lead_variant_id", "chr", "pos_hg38", "ref", "alt", "peak_id",
          "peak_start_hg38", "peak_stop_hg38", "heldout_fold"]
LENGTHS = (2048, 16384, 1048576)
PILOT_SIZE = 64
PILOT_SEED = 20260930
FASTA = FAI.with_suffix("")
CHECKPOINT = ROOT / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-weights-20260915T103442Z/checkpoints"
SCORING_PYTHON = Path("/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/bin/python")
PILOT_RULE = "Eligible whole-target-contained 1Mb rows only; strata are chromosome fold crossed with zero-mask-overlap, positive-partial-peak-overlap or whole-peak-inside-mask. Sort each stratum by SHA256 of seed TAB canonical_variant_key TAB peak_id; cycle sorted stratum names taking one unused row per stratum until64 rows or exhaustion. No beta/p-value/sequence prediction is used."
RESEARCH_SOURCE = ROOT / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-runtime-probe-20260915T104500Z/code/alphagenome_research/src/alphagenome_research/model"
SDK_SOURCE = Path("/gpfs/commons/home/jameslee/micromamba/envs/alphagenome_local/lib/python3.11/site-packages/alphagenome/data")


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def chrom(value):
    return str(value).removeprefix("chr")


def integer(value):
    try:
        x = float(value)
        return int(x) if np.isfinite(x) and x == int(x) else None
    except (TypeError, ValueError, OverflowError):
        return None


def overlap(a, b, c, d):
    return max(0, min(b, d) - max(a, c))


def write_json(path, value):
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write("\n")


def main():
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run this census in its authorized SLURM job.")
    output = Path(os.environ["CODEX_REC_OUTPUT"]) / ("target_region_census_" + os.environ["SLURM_JOB_ID"])
    output.mkdir(parents=True, exist_ok=False)
    inputs = {"labels": LABELS, "peaks": PEAKS, "reference_fai": FAI,
              "script": Path(__file__), "launcher": Path(__file__).with_name("run_target_region_census.sbatch"),
              "installed_native_dna_model": RESEARCH_SOURCE / "dna_model.py",
              "installed_center_mask": RESEARCH_SOURCE / "variant_scoring/center_mask.py",
              "installed_track_data": SDK_SOURCE / "track_data.py",
              "installed_genome_coordinates": SDK_SOURCE / "genome.py"}
    receipt = {
        "inputs": {name: {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
                   for name, path in inputs.items()},
        "python": sys.version, "platform": platform.platform(),
        "numpy": np.__version__, "pandas": pd.__version__,
        "metadata_columns_read": FIELDS,
        "peak_columns_read": ["#chr", "start", "end", "peakID"],
        "folds": [1, 2, 3, 4], "outcome_columns_loaded": False,
        "coordinate_semantics": "GRCh38 BED0 half-open peaks; one-based variant; variant input window start=position1-1-length//2",
        "assembly_check_limit": "Independent deposited peak coordinates compared to explicit hg38 C2 columns; no sequence/genotype verification or per-row assembly label supplied.",
        "native_output_semantics": "Installed predict_sequence passes the provided genomic interval unchanged to TrackData; TrackData validates resolution*positional_bins equals interval width. At 1bp, successful REF/ALT output therefore covers the entire supplied interval with no hidden coordinate crop.",
        "center_mask_semantics": "Installed create_center_mask at resolution1,width501 uses [variant.start-250,variant.start+251), clipped to supplied interval; Variant.start=position1-1. All assessed window lengths exceed501 and place the variant at length//2.",
        "model_output_limit": "Source inspection establishes output coordinate semantics conditional on successful inference; alphabet/REF validation and actual returned shape/interval/track identities still require scoring-time checks. Containment does not establish accuracy at model edges.",
        "seed": None, "stochastic_operations": False,
        "pilot_selection": {"size": PILOT_SIZE, "hash_seed": PILOT_SEED, "rule": PILOT_RULE},
        "native_resources": {
            "fasta": {"path": str(FASTA), "exists": FASTA.is_file(),
                      "bytes": FASTA.stat().st_size if FASTA.is_file() else None,
                      "hash_state": "FAI hashed; full FASTA hash not computed by this metadata-only census"},
            "checkpoint": {"path": str(CHECKPOINT), "exists": CHECKPOINT.is_dir()},
            "scoring_python": {"path": str(SCORING_PYTHON), "exists": SCORING_PYTHON.is_file()},
        },
    }
    write_json(output / "receipt.json", receipt)
    # usecols prevents loading every effect, p-value or phenotype column,
    # including those belonging to closed fold 0.
    labels = pd.read_csv(LABELS, sep="\t", usecols=FIELDS, dtype=str, keep_default_na=False)
    folds = pd.to_numeric(labels.heldout_fold, errors="raise").astype(int)
    if not folds.isin(range(5)).all():
        raise ValueError("Unexpected chromosome fold.")
    metadata_rows_total = len(labels)
    labels = labels.loc[folds.isin([1, 2, 3, 4])].copy()
    labels["heldout_fold"] = folds.loc[labels.index]
    if labels.empty or labels.lead_variant_id.duplicated().any():
        raise ValueError("C2 development variant identities missing or duplicated.")
    labels["chromosome"] = labels["chr"].map(chrom)
    if labels.groupby("chromosome").heldout_fold.nunique().max() != 1:
        raise ValueError("A chromosome crosses development folds.")
    peaks = pd.read_csv(PEAKS, sep="\t", usecols=["#chr", "start", "end", "peakID"],
                        dtype=str, keep_default_na=False)
    peaks["source_chromosome"] = peaks["#chr"].map(chrom)
    geometries = peaks[["peakID", "source_chromosome", "start", "end"]].drop_duplicates()
    if geometries.peakID.duplicated().any():
        write_json(output / "failure.json", {"reason": "ambiguous_duplicate_peak_geometry",
                   "ambiguous_ids": int(geometries.loc[geometries.peakID.duplicated(False), "peakID"].nunique())})
        raise ValueError("One source peak ID maps to conflicting geometry; no nearest-peak repair allowed.")
    if (geometries.peakID == "").any():
        raise ValueError("Blank source peak identifier.")
    axis = geometries.set_index("peakID").to_dict("index")
    lengths = {}
    reference_aliases = {}
    with FAI.open() as handle:
        for line in handle:
            parts = line.rstrip().split("\t")
            key = chrom(parts[0])
            if key in lengths and lengths[key] != int(parts[1]):
                raise ValueError("Reference chromosome alias has conflicting length.")
            if key in reference_aliases and reference_aliases[key] != parts[0]:
                raise ValueError("Ambiguous reference chromosome alias; do not choose arbitrarily.")
            lengths[key] = int(parts[1])
            reference_aliases[key] = parts[0]
    # Half-open edge checks are scientific invariants, evaluated in the job.
    assert overlap(0, 10, 10, 20) == 0
    assert overlap(0, 10, 9, 20) == 1
    assert overlap(2, 8, 0, 10) == 6
    rows = []
    for r in labels.itertuples(index=False):
        pos = integer(r.pos_hg38)
        start, end = integer(r.peak_start_hg38), integer(r.peak_stop_hg38)
        source = axis.get(r.peak_id)
        reason = "geometry_agrees"
        source_start = integer(source["start"]) if source else None
        source_end = integer(source["end"]) if source else None
        parts = r.lead_variant_id.split(":")
        if len(parts) != 4 or chrom(parts[0]) != r.chromosome or integer(parts[1]) != pos or parts[2].upper() != r.ref.upper() or parts[3].upper() != r.alt.upper():
            reason = "variant_identifier_mismatch"
        elif source is None:
            reason = "unmapped_peak_id"
        elif None in (pos, start, end, source_start, source_end):
            reason = "invalid_coordinate"
        elif r.chromosome not in lengths:
            reason = "unknown_reference_chromosome"
        elif source["source_chromosome"] != r.chromosome:
            reason = "target_chromosome_mismatch"
        elif (start, end) != (source_start, source_end):
            reason = "source_hg38_coordinate_disagreement"
        elif not (0 <= start < end <= lengths[r.chromosome]) or not (1 <= pos <= lengths[r.chromosome]):
            reason = "invalid_reference_bounds"
        valid = reason == "geometry_agrees"
        for length in LENGTHS:
            ws = pos - 1 - length // 2 if pos is not None else None
            we = ws + length if ws is not None else None
            window_valid = bool(valid and ws >= 0 and we <= lengths[r.chromosome])
            contained = bool(window_valid and ws <= start and end <= we)
            state = reason if not valid else ("window_out_of_reference_bounds" if not window_valid else
                    ("whole_peak_contained" if contained else ("partial_peak_in_window" if overlap(ws, we, start, end) else "peak_outside_window")))
            mask_start = pos - 1 - 250 if valid else None
            mask_end = mask_start + 501 if valid else None
            if valid:
                assert mask_end - mask_start == 501
                assert ws <= mask_start < mask_end <= we
            mask_overlap = overlap(mask_start, mask_end, start, end) if valid else None
            rows.append({
                "variant_id": r.lead_variant_id, "peak_id": r.peak_id, "chromosome": r.chromosome,
                "heldout_fold": int(r.heldout_fold), "length_bp": length,
                "peak_start0": start, "peak_end0": end,
                "source_peak_start0": source_start, "source_peak_end0": source_end,
                "coordinate_state": reason, "target_readout_state": state,
                "window_start0": ws, "window_end0": we, "whole_target_in_input": contained,
                "peak_width_bp": end - start if valid else None,
                "variant_inside_peak": bool(start <= pos - 1 < end) if valid else None,
                "signed_variant_to_peak_center_bp": pos - 1 - (start + end) / 2 if valid else None,
                "mask_start0": mask_start, "mask_end0": mask_end, "mask_overlap_bp": mask_overlap,
                "mask_overlap_fraction": mask_overlap / 501 if valid else None,
                "peak_overlap_fraction": mask_overlap / (end - start) if valid else None,
            })
    table = pd.DataFrame(rows)
    table.to_csv(output / "development_target_geometry.tsv.gz", sep="\t", index=False, mode="x")
    counts = table.groupby(["heldout_fold", "length_bp", "coordinate_state", "target_readout_state"], dropna=False).size().reset_index(name="rows")
    counts.to_csv(output / "development_geometry_counts.tsv", sep="\t", index=False, mode="x")
    overlap_summaries = []
    for (fold, length), group in table.groupby(["heldout_fold", "length_bp"]):
        for scope, subset in (
            ("geometry_agrees", group.loc[group.coordinate_state == "geometry_agrees"]),
            ("whole_peak_and_window_admitted", group.loc[group.target_readout_state == "whole_peak_contained"]),
        ):
            overlap_summaries.append({
                "heldout_fold": int(fold), "length_bp": int(length), "scope": scope,
                "rows": len(subset),
                "mask_zero_overlap": int((subset.mask_overlap_bp == 0).sum()),
                "mask_positive_overlap": int((subset.mask_overlap_bp > 0).sum()),
                "whole_peak_inside_501_mask": int((subset.mask_overlap_bp == subset.peak_width_bp).sum()),
                "whole_501_mask_inside_peak": int((subset.mask_overlap_bp == 501).sum()),
                "variant_inside_peak": int((subset.variant_inside_peak == True).sum()),
            })
    pd.DataFrame(overlap_summaries).to_csv(output / "development_mask_overlap_counts.tsv", sep="\t", index=False, mode="x")
    # Produce inference inputs with no beta/p-value columns. Chromosome aliases
    # are exact contig names from the same FAI used by the current native scorer.
    admitted = table.loc[(table.length_bp == 1048576) &
                         (table.target_readout_state == "whole_peak_contained")].copy()
    manifest = admitted.merge(labels[["lead_variant_id", "pos_hg38", "ref", "alt"]],
                              left_on="variant_id", right_on="lead_variant_id", validate="one_to_one")
    manifest["key"] = manifest.variant_id.str.removeprefix("chr")
    manifest["chr"] = manifest.chromosome.map(reference_aliases)
    manifest["assembly"] = "GRCh38"
    manifest["reference_fai_sha256"] = receipt["inputs"]["reference_fai"]["sha256"]
    manifest["source_labels_sha256"] = receipt["inputs"]["labels"]["sha256"]
    manifest["source_peaks_sha256"] = receipt["inputs"]["peaks"]["sha256"]
    manifest["readout_population"] = "historically_significance_selected_Currin_variant_peak_associations;not_gene_effects_or_arbitrary_target_predictions"
    manifest["sequence_validation_state"] = "pending_scoring_time_REF_alphabet_and_runtime_checks"
    manifest["pilot_geometry"] = np.where(manifest.mask_overlap_bp == 0, "zero_mask_overlap",
        np.where(manifest.mask_overlap_bp == manifest.peak_width_bp, "whole_peak_inside_mask", "positive_partial_peak_overlap"))
    manifest["pilot_stratum"] = manifest.heldout_fold.astype(str) + ":" + manifest.pilot_geometry
    manifest["pilot_hash"] = [hashlib.sha256(f"{PILOT_SEED}\t{key}\t{peak}".encode()).hexdigest()
                               for key, peak in zip(manifest.key, manifest.peak_id)]
    columns = ["key", "variant_id", "chr", "pos_hg38", "ref", "alt", "peak_id", "assembly", "heldout_fold",
               "length_bp", "window_start0", "window_end0", "peak_start0", "peak_end0", "mask_start0", "mask_end0",
               "mask_overlap_bp", "peak_width_bp", "pilot_geometry", "pilot_stratum", "pilot_hash",
               "coordinate_state", "target_readout_state", "sequence_validation_state", "readout_population",
               "reference_fai_sha256", "source_labels_sha256", "source_peaks_sha256"]
    manifest = manifest[columns].sort_values(["heldout_fold", "chr", "pos_hg38", "peak_id"]).reset_index(drop=True)
    if manifest.empty or manifest.key.duplicated().any() or manifest["chr"].isna().any():
        raise ValueError("Invalid admitted target manifest.")
    # Round-robin hash selection balances geometry and folds without effects.
    queues = {name: list(group.sort_values(["pilot_hash", "key", "peak_id"]).index)
              for name, group in manifest.groupby("pilot_stratum")}
    selected = []
    depth = 0
    while len(selected) < min(PILOT_SIZE, len(manifest)):
        progressed = False
        for name in sorted(queues):
            if depth < len(queues[name]):
                selected.append(queues[name][depth])
                progressed = True
                if len(selected) == min(PILOT_SIZE, len(manifest)):
                    break
        if not progressed:
            break
        depth += 1
    pilot = manifest.loc[selected].copy()
    pilot.insert(0, "pilot_selection_order", range(len(pilot)))
    assert len(pilot) == min(PILOT_SIZE, len(manifest)) and not pilot.key.duplicated().any()
    manifest_path = output / "development_target_manifest_1mb.tsv.gz"
    pilot_path = output / "target_readout_pilot64.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False, mode="x")
    pilot.to_csv(pilot_path, sep="\t", index=False, mode="x")
    rejected = table.loc[(table.length_bp == 1048576) & (table.target_readout_state != "whole_peak_contained")]
    rejected[["variant_id", "peak_id", "heldout_fold", "coordinate_state", "target_readout_state"]].to_csv(
        output / "target_manifest_geometry_exclusions.tsv", sep="\t", index=False, mode="x")
    first = table.loc[table.length_bp == LENGTHS[0]]
    summary = {
        "status": "metadata_census_complete", "metadata_rows_total": metadata_rows_total,
        "development_rows": len(labels), "closed_fold_metadata_rows_excluded": metadata_rows_total - len(labels),
        "source_peak_rows": len(peaks), "source_unique_geometries": len(geometries),
        "identical_source_peak_rows_collapsed": len(peaks) - len(geometries),
        "development_distinct_targets": int(labels.peak_id.nunique()),
        "repeated_development_targets": int((labels.peak_id.value_counts() > 1).sum()),
        "variant_target_pairs": len(labels[["lead_variant_id", "peak_id"]].drop_duplicates()),
        "coordinate_states": {str(k): int(v) for k, v in first.coordinate_state.value_counts().items()},
        "mask_overlap_counts_by_fold_length": overlap_summaries,
        "development_fold_rows": {str(k): int(v) for k, v in labels.heldout_fold.value_counts().sort_index().items()},
        "target_manifest": {"path": str(manifest_path), "rows": len(manifest), "sha256": sha256(manifest_path),
                            "geometry_exclusions": len(rejected)},
        "pilot_manifest": {"path": str(pilot_path), "rows": len(pilot), "sha256": sha256(pilot_path),
                           "hash_seed": PILOT_SEED, "selection_rule": PILOT_RULE,
                           "strata": {str(k): int(v) for k, v in pilot.pilot_stratum.value_counts().sort_index().items()}},
        "assembly_mismatch_count": None,
        "assembly_mismatch_reason": "No per-row assembly annotation; coordinate disagreements counted, without inferring their cause.",
        "lengths": {str(length): {str(k): int(v) for k, v in table.loc[table.length_bp == length, "target_readout_state"].value_counts().items()} for length in LENGTHS},
        "half_open_boundary_checks": "passed", "outcomes_read": False,
        "limit": receipt["model_output_limit"],
    }
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
