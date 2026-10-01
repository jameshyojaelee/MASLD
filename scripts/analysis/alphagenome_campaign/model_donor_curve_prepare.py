#!/usr/bin/env python3
"""Outcome-independent nested counted-interval panels and paired donor inputs."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time

import numpy as np
import pandas as pd
import pysam

import review_donor as D
PROJ = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJ/"scripts/analysis/alphagenome_program"))
import i1_common as C

BASE = PROJ/"GWAS/finemapping/results/alphagenome_campaign/week1-20260915/model"
ORIGINAL = BASE/"donor_fixture_21773589/regions.tsv"
SIZES = (128, 1024, 8192)
LENGTH = 2048
SEED = 20260915


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")


def hash_key(value):
    return hashlib.sha256(str(value).encode()).hexdigest()


def window_components(frame):
    """Connected components of actual input intervals, including transitivity."""
    output = pd.Series(index=frame.index, dtype=str)
    for chrom, group in frame.groupby("chrom", sort=True):
        end, component = -1, None
        for index, row in group.sort_values(["window_start0", "window_end0", "region_key"]).iterrows():
            if row.window_start0 >= end:
                component = f"{chrom}:{int(row.window_start0)}"
                end = int(row.window_end0)
            else:
                end = max(end, int(row.window_end0))
            output.loc[index] = component
    return output


def make_manifest(args):
    start = time.monotonic()
    args.out.mkdir(parents=True, exist_ok=False)
    axis_path = D.FIX/"molecular/h3k27ac_feature_axis.tsv"
    source = pd.read_csv(axis_path, sep="\t")
    if len(source) != 96460 or source.opaque_source_feature_key.duplicated().any():
        raise ValueError("Full H3 region axis differs")
    if not np.array_equal(source.h3k27ac_feature_index, np.arange(96460)):
        raise ValueError("H3 region indices are not the frozen count axis")
    original = pd.read_csv(args.original, sep="\t")
    if len(original) != 128 or original.groupby("region_role").size().to_dict() != {"held":64, "train":64}:
        raise ValueError("Original pilot definition changed")
    records = []
    fasta = pysam.FastaFile(C.FASTA_PATH)
    chromosome_lengths = dict(zip(fasta.references, fasta.lengths))
    for item in source.itertuples():
        key = str(item.opaque_source_feature_key)
        match = re.fullmatch(r"(chr[^:]+):(\d+)-(\d+)", key)
        if match is None:
            raise ValueError(f"Unparseable counted interval identity: {key}")
        chrom, left, right = match.groups()
        lo, hi = int(left)-1, int(right)
        if lo < 0 or hi <= lo:
            raise ValueError("Invalid source counted interval")
        center = (lo+hi)//2
        a, b = center-LENGTH//2, center+LENGTH//2
        reason = "eligible"
        sequence_hash = ""
        if chrom not in {f"chr{i}" for i in range(1,23)}:
            reason = "nonautosomal"
        elif hi-lo > LENGTH:
            reason = "counted_interval_wider_than2048"
        elif chrom not in chromosome_lengths or a < 0 or b > chromosome_lengths[chrom]:
            reason = "incomplete_reference_window"
        elif lo < a or hi > b:
            reason = "whole_counted_interval_not_contained"
        else:
            sequence = fasta.fetch(chrom, a, b).upper()
            if len(sequence) != LENGTH or set(sequence)-set("ACGT"):
                reason = "nonACGT_or_incomplete_window"
            else:
                sequence_hash = hash_key(sequence)
        records.append({"region_index": int(item.h3k27ac_feature_index), "region_key": key,
            "chrom": chrom, "start0": lo, "end0": hi, "counted_width": hi-lo,
            "window_start0": a, "window_end0": b, "eligibility": reason,
            "sequence_sha256": sequence_hash, "coordinate_residual_bp": 1,
            "coordinate_semantics": "counted_interval", "sample_hash": hash_key("b2-region-curve|"+key)})
    population = pd.DataFrame(records)
    eligible = population[population.eligibility == "eligible"].copy()
    eligible["input_component"] = window_components(eligible)
    old_roles = original.set_index("region_key").region_role.to_dict()
    eligible["original_role"] = eligible.region_key.map(old_roles).fillna("not_original")
    old_held = eligible[eligible.original_role == "held"]
    reserved = set(old_held.input_component)
    old_train_components = set(eligible.loc[eligible.original_role == "train", "input_component"])-reserved
    held_keys = list(old_held.region_key)
    # Protect compatible old training anchors while fixing the new held panel.
    # Whole components containing held rows are reserved from every train size.
    for row in eligible.sort_values("sample_hash").itertuples():
        if len(held_keys) == args.held_regions:
            break
        if row.region_key in held_keys or row.input_component in old_train_components:
            continue
        held_keys.append(row.region_key)
        reserved.add(row.input_component)
    if len(held_keys) != args.held_regions:
        raise ValueError("Insufficient eligible held regions under component separation")
    training = eligible[~eligible.input_component.isin(reserved)].copy()
    training["anchor_priority"] = (training.original_role != "train").astype(int)
    training = training.sort_values(["anchor_priority", "sample_hash"]).head(max(SIZES)).copy()
    if len(training) != max(SIZES):
        raise ValueError("Insufficient training regions after complete component holdout")
    training["training_rank"] = np.arange(1, len(training)+1)
    training["region_role"] = "train"
    held = eligible.set_index("region_key").loc[held_keys].reset_index().copy()
    held["training_rank"] = -1
    held["region_role"] = "held"
    chosen = pd.concat([training, held], ignore_index=True)
    chosen["inner_region_role"] = ["outer_held" if role == "held" else (
        "inner_valid" if int(hash_key("b2-curve-inner|"+component)[:16],16)%5 == 0 else "inner_train")
        for role, component in zip(chosen.region_role, chosen.input_component)]
    chosen["fixed_trained_evaluation"] = (chosen.region_role == "train") & (chosen.training_rank <=128)
    for size in SIZES:
        chosen[f"training_{size}"] = (chosen.region_role == "train") & (chosen.training_rank <= size)
        selected = chosen[chosen[f"training_{size}"]]
        if len(selected) != size or set(selected.inner_region_role) != {"inner_train", "inner_valid"}:
            raise ValueError("Nested size lacks both inner-region partitions")
        if selected.groupby("inner_region_role").input_component.nunique().min() < 4:
            raise ValueError("Fewer than four inner independent input components")
    if chosen.groupby("input_component").region_role.nunique().max() != 1:
        raise ValueError("An actual input-overlap component crosses outer region holdout")
    if chosen.groupby("input_component").inner_region_role.nunique().max() != 1:
        raise ValueError("An actual input-overlap component crosses inner region holdout")
    # All old region identities must remain tied to the original integer axis.
    check = original.merge(population, on="region_key", suffixes=("_old", "_new"), validate="one_to_one")
    for column in ("region_index", "chrom", "start0", "end0"):
        if not np.array_equal(check[column+"_old"], check[column+"_new"]):
            raise ValueError("Old counted-interval identity changed")
    population.to_csv(args.out/"full_region_eligibility.tsv.gz", sep="\t", index=False)
    chosen.to_csv(args.out/"regions.tsv", sep="\t", index=False)
    dispositions = original[["region_key", "region_role"]].merge(
        chosen[["region_key", "region_role", "input_component"]], on="region_key", how="left", suffixes=("_old", "_new"))
    dispositions["disposition"] = np.where(dispositions.region_role_new.notna(), "retained", "excluded_or_component_conflict")
    dispositions.to_csv(args.out/"original128_disposition.tsv", sep="\t", index=False)
    counts = population.groupby("eligibility").size().to_dict()
    receipt = {"source_regions": len(population), "eligible_regions": len(eligible), "eligibility_counts": counts,
        "selected_total": len(chosen), "training_sizes": SIZES, "fixed_held_regions": len(held),
        "original_held_retained": int((held.original_role == "held").sum()),
        "original_training_anchors_retained": int((training.original_role == "train").sum()),
        "input_overlap_components": int(chosen.input_component.nunique()),
        "inner_counts": {str(size): chosen[chosen[f"training_{size}"]].inner_region_role.value_counts().to_dict() for size in SIZES},
        "input_length": LENGTH, "coordinate_residual_bp": 1, "coordinate_semantics": "counted_interval",
        "outcomes_read": False, "count_matrix_loaded": False, "fasta": C.FASTA_PATH,
        "axis_sha256": D.sha(axis_path), "manifest_sha256": D.sha(args.out/"regions.tsv"),
        "original_manifest_sha256": D.sha(args.original), "code_sha256": D.sha(__file__),
        "selection": "Coordinate hash after full input eligibility; complete actual input-overlap components separate outer and inner held regions",
        "seconds": time.monotonic()-start}
    write(args.out/"manifest_receipt.json", receipt)
    print(json.dumps(receipt), flush=True)


def paired(args):
    reg = pd.read_csv(args.out/"regions.tsv", sep="\t")
    receipt = json.loads((args.out/"manifest_receipt.json").read_text())
    if receipt["outcomes_read"] or receipt["manifest_sha256"] != D.sha(args.out/"regions.tsv"):
        raise ValueError("Outcome-free manifest definition differs")
    if (args.out/"paired.npz").exists():
        raise FileExistsError(args.out/"paired.npz")
    ids = pd.read_csv(D.FIX/"molecular/participant_axis.tsv", sep="\t").participant_id.astype(str).to_numpy()
    ft = pd.read_csv(D.FIX/"folds/participant_outer_folds.tsv", sep="\t").set_index("participant_id")
    fold = ft.loc[ids,"outer_fold"].to_numpy(int)
    raw_x = np.load(D.FIX/"molecular/rna_values.npy", mmap_mode="r")
    raw_y = np.load(D.FIX/"molecular/h3k27ac_counts.npy", mmap_mode="r")
    if len(ids) != 99 or len(set(ids)) !=99 or raw_x.shape != (99,42163) or raw_y.shape != (99,96460):
        raise ValueError("Paired molecular source identities/dimensions changed")
    if not np.issubdtype(raw_y.dtype, np.integer):
        raise ValueError("Source H3 counts lost integer units")
    for name in ("rna", "h3k27ac"):
        if not np.load(D.FIX/f"molecular/{name}_observed_mask.npy").all():
            raise ValueError("Missing source assays cannot become zero targets")
    x = D.logcpm(raw_x)
    y = D.logcpm(raw_y)[:,reg.region_index.to_numpy(int)]
    genes = pd.read_csv(D.FIX/"molecular/rna_feature_axis.tsv", sep="\t").stable_gene_id.to_numpy(str)
    np.savez_compressed(args.out/"paired.npz", participant_ids=ids, donor_fold=fold,
        rna_log2cpm_full_library=x, target_log2cpm_full_library=y, gene_ids=genes,
        region_key=reg.region_key.to_numpy(str), full_H3_library_totals=raw_y.sum(1,dtype=np.float64),
        full_RNA_library_totals=raw_x.sum(1,dtype=np.float64))
    write(args.out/"paired_receipt.json", {"participants":99,"source_regions":96460,"selected_regions":len(reg),
        "RNA_features":42163,"normalization":"log2(1+CPM), full source feature library before subsetting",
        "continuous_RNA_not_integer_counts":True,"no_outcome_selected_regions":True,
        "manifest_sha256":D.sha(args.out/"regions.tsv"),"paired_sha256":D.sha(args.out/"paired.npz"),
        "assay_pairing":"same participant, same sample, different aliquots", "source_data_rehosting":False})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode",choices=("manifest","paired"))
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--original",type=Path,default=ORIGINAL)
    parser.add_argument("--held-regions",type=int,default=1024)
    arguments = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute allocation required")
    (make_manifest if arguments.mode == "manifest" else paired)(arguments)
