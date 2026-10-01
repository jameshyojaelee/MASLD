"""Build an outcome-free 1Mb target manifest from a completed geometry census."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
FIELDS = ["lead_variant_id", "chr", "pos_hg38", "ref", "alt", "peak_id", "heldout_fold"]
PILOT_SIZE = 64
PILOT_SEED = 20260930
RULE = "Whole-target-contained 1Mb development folds1--4; fold crossed with zero501-mask overlap, positive overlap without full peak, or full peak inside501mask; SHA256(seed TAB canonical_variant_key TAB peak_id) order within each stratum; round-robin sorted stratum names until64 rows. No outcomes or predictions."


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main(census):
    if not os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Requires authorized compute job.")
    out = Path(os.environ["CODEX_REC_OUTPUT"]) / ("target_region_manifest_" + os.environ["SLURM_JOB_ID"])
    out.mkdir(parents=True, exist_ok=False)
    receipt = json.loads((census / "receipt.json").read_text())
    summary = json.loads((census / "summary.json").read_text())
    if summary["status"] != "metadata_census_complete" or summary["outcomes_read"]:
        raise ValueError("Source census is not a completed outcome-free result.")
    for name in ("labels", "peaks", "reference_fai"):
        if sha256(receipt["inputs"][name]["path"]) != receipt["inputs"][name]["sha256"]:
            raise ValueError("Census input changed: " + name)
    design = {"source_census": str(census), "source_receipt_sha256": sha256(census / "receipt.json"),
              "source_summary_sha256": sha256(census / "summary.json"), "seed": PILOT_SEED,
              "pilot_size": PILOT_SIZE, "selection_rule": RULE, "metadata_columns": FIELDS,
              "outcomes_read": False, "python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
              "script_sha256": sha256(__file__),
              "sequence_validation": "REF/alphabet/runtime validation pending; no model imported or run."}
    # Record the complete selection rule before constructing a pilot.
    (out / "design.json").write_text(json.dumps(design, indent=2) + "\n")
    geometry_cols = ["variant_id", "peak_id", "chromosome", "heldout_fold", "length_bp", "coordinate_state",
                     "target_readout_state", "window_start0", "window_end0", "peak_start0", "peak_end0",
                     "mask_start0", "mask_end0", "mask_overlap_bp", "peak_width_bp"]
    geometry_path = census / "development_target_geometry.tsv.gz"
    geometry = pd.read_csv(geometry_path, sep="\t", usecols=geometry_cols)
    geometry = geometry.loc[geometry.length_bp == 1048576].copy()
    if len(geometry) != summary["development_rows"] or not geometry.heldout_fold.isin([1, 2, 3, 4]).all():
        raise ValueError("Source geometry population/folds disagree.")
    if geometry.variant_id.duplicated().any():
        raise ValueError("Duplicate geometry variant; arbitrary target repair prohibited.")
    geometry.loc[geometry.target_readout_state != "whole_peak_contained"].to_csv(
        out / "geometry_exclusions.tsv", sep="\t", index=False, mode="x")
    admitted = geometry.loc[geometry.target_readout_state == "whole_peak_contained"].copy()
    if not admitted.coordinate_state.eq("geometry_agrees").all():
        raise ValueError("Admitted target does not independently match source geometry.")
    labels = pd.read_csv(receipt["inputs"]["labels"]["path"], sep="\t", usecols=FIELDS, dtype=str)
    labels = labels.loc[labels.heldout_fold.isin(["1", "2", "3", "4"])].copy()
    aliases = {}
    with Path(receipt["inputs"]["reference_fai"]["path"]).open() as handle:
        for line in handle:
            name = line.split("\t")[0]
            key = name.removeprefix("chr")
            if key in aliases and aliases[key] != name:
                raise ValueError("Ambiguous reference contig alias.")
            aliases[key] = name
    manifest = admitted.merge(labels[["lead_variant_id", "pos_hg38", "ref", "alt"]],
                               left_on="variant_id", right_on="lead_variant_id", validate="one_to_one")
    if len(manifest) != len(admitted):
        raise ValueError("Metadata join lost eligible rows.")
    manifest["key"] = manifest.variant_id.str.removeprefix("chr")
    manifest["chr"] = manifest.chromosome.astype(str).map(aliases)
    manifest["assembly"] = "GRCh38"
    manifest["sequence_validation_state"] = "pending_REF_alphabet_runtime_checks"
    manifest["population"] = "historically_significance_selected_Currin_variant_peak_associations;not_gene_effects_or_arbitrary_targets"
    manifest["geometry_stratum"] = np.where(manifest.mask_overlap_bp == 0, "zero_mask_overlap",
        np.where(manifest.mask_overlap_bp == manifest.peak_width_bp, "whole_peak_inside_mask", "positive_partial_peak_overlap"))
    manifest["pilot_stratum"] = manifest.heldout_fold.astype(str) + ":" + manifest.geometry_stratum
    manifest["pilot_hash"] = [hashlib.sha256(f"{PILOT_SEED}\t{key}\t{peak}".encode()).hexdigest()
                               for key, peak in zip(manifest.key, manifest.peak_id)]
    for name in ("labels", "peaks", "reference_fai"):
        manifest[name + "_sha256"] = receipt["inputs"][name]["sha256"]
    manifest = manifest.drop(columns="lead_variant_id").sort_values(["heldout_fold", "chr", "pos_hg38", "peak_id"]).reset_index(drop=True)
    if manifest.empty or manifest.key.duplicated().any() or manifest["chr"].isna().any():
        raise ValueError("Invalid outcome-free target manifest.")
    queues = {name: list(group.sort_values(["pilot_hash", "key", "peak_id"]).index)
              for name, group in manifest.groupby("pilot_stratum")}
    picked = []
    depth = 0
    while len(picked) < min(PILOT_SIZE, len(manifest)):
        for name in sorted(queues):
            if depth < len(queues[name]):
                picked.append(queues[name][depth])
                if len(picked) == min(PILOT_SIZE, len(manifest)):
                    break
        depth += 1
    pilot = manifest.loc[picked].copy()
    pilot.insert(0, "pilot_selection_order", range(len(pilot)))
    manifest_path, pilot_path = out / "development_target_manifest_1mb.tsv.gz", out / "target_readout_pilot64.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False, mode="x")
    pilot.to_csv(pilot_path, sep="\t", index=False, mode="x")
    result = {"status": "outcome_free_manifest_complete", "rows": len(manifest), "pilot_rows": len(pilot),
              "source_geometry_sha256": sha256(geometry_path), "manifest_sha256": sha256(manifest_path),
              "pilot_sha256": sha256(pilot_path), "geometry_exclusions": len(geometry) - len(manifest),
              "pilot_strata": {str(k): int(v) for k, v in pilot.pilot_stratum.value_counts().sort_index().items()},
              "outcomes_read": False, "sequence_inference_run": False}
    (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", type=Path, required=True)
    main(parser.parse_args().census)
