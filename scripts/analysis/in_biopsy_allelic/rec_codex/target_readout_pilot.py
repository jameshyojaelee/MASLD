"""Outcome-free same-prediction center versus associated-peak ATAC readouts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "scripts/analysis/alphagenome_program"))
import i1_common as C

PILOT_SHA = "577901e5058d0238d3837fcadd8497ec58b610972469363fea70fe4d954180e8"
FULL_SHA = "ef47cc5b60aa2aab542c686772a620575af3f90e1567eac7c17beb1ac2c63eef"
FULL_ROWS = 25477
FIELDS = ["pilot_selection_order", "key", "chr", "pos_hg38", "ref", "alt", "heldout_fold",
          "peak_id", "peak_start0", "peak_end0", "window_start0", "window_end0",
          "mask_start0", "mask_end0", "geometry_stratum"]


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def peak_mask(interval, start, end):
    if not interval.start <= start < end <= interval.end:
        raise C.ContractError("Whole associated peak is not contained; clipping prohibited.")
    mask = np.zeros((interval.width, 1), dtype=bool)
    mask[start-interval.start:end-interval.start] = True
    return mask


def read_manifest(manifest, full_development=False):
    expected_sha = FULL_SHA if full_development else PILOT_SHA
    expected_rows = FULL_ROWS if full_development else 64
    fields = FIELDS[1:] if full_development else FIELDS
    if C.sha256_file(manifest) != expected_sha:
        raise C.ContractError("Manifest differs from its fixed metadata-only population.")
    rows = pd.read_csv(manifest, sep="\t", usecols=fields)
    order = "key" if full_development else "pilot_selection_order"
    rows = rows.sort_values(order)
    if len(rows) != expected_rows or rows.key.duplicated().any() or not rows.heldout_fold.isin([1, 2, 3, 4]).all():
        raise C.ContractError("Fixed population/folds invalid.")
    return rows, expected_sha


def main(manifest, full_development=False):
    if not os.environ.get("SLURM_JOB_ID"):
        raise C.ContractError("SLURM compute allocation required.")
    rows, manifest_sha = read_manifest(manifest, full_development)
    prefix = "target_readout_development_" if full_development else "target_readout_pilot_"
    out = Path(os.environ["CODEX_REC_OUTPUT"]) / (prefix + os.environ["SLURM_JOB_ID"])
    out.mkdir(parents=True, exist_ok=False)
    card = C.gpu_card()
    if card.get("card_tag") != "l40s":
        raise C.ContractError("Requires the existing L40S native scoring route.")
    import jax
    import jax.numpy as jnp
    import pysam
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    scorer, _ = C.atac_dnase_scorers()
    make_mask, aggregate, _, resolution = C.center_mask_and_aggregation()
    if resolution(scorer.requested_output) != 1:
        raise C.ContractError("ATAC must have 1-bp resolution.")
    # Small numerical invariants on the allocated device, before model restoration.
    test_interval = genome.Interval(chromosome="chr1", start=100, end=108)
    test_mask = peak_mask(test_interval, 102, 106)
    if np.flatnonzero(test_mask[:, 0]).tolist() != [2, 3, 4, 5]:
        raise C.ContractError("BED0 half-open mask failed.")
    for bad_start, bad_end in [(99, 106), (102, 109), (106, 106)]:
        try:
            peak_mask(test_interval, bad_start, bad_end)
        except C.ContractError:
            pass
        else:
            raise C.ContractError("Partial/empty mask accepted.")
    x = np.arange(24, dtype=np.float32).reshape(8, 3)
    y = x + np.array([0, 2, 5], dtype=np.float32)
    got = np.asarray(aggregate(jnp.asarray(x), jnp.asarray(y), jnp.asarray(test_mask),
                              aggregation_type=scorer.aggregation_type))
    expected = np.log2(1+y[2:6].sum(axis=0))-np.log2(1+x[2:6].sum(axis=0))
    if not np.allclose(got, expected, atol=2e-6, rtol=0):
        raise C.ContractError("Installed aggregation differs from independent sum-then-log formula.")
    receipt = {"manifest": str(manifest), "pilot_sha256": PILOT_SHA,
               "manifest_sha256": manifest_sha, "selected_rows": len(rows),
               "run_mode": "full_development" if full_development else "pilot",
               "script_sha256": C.sha256_file(__file__), "helper_sha256": C.sha256_file(C.__file__),
               "card": card, "python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
               "jax": jax.__version__, "checkpoint": str(C.CHECKPOINT), "fasta": C.FASTA_PATH,
               "outcomes_read": False, "requested_outputs": ["ATAC", "DNASE"],
               "tracks": C.AG_ATAC_TRACKS, "aggregation": "DIFF_LOG2_SUM",
               "mask_rule": "501bp variant center and complete associated BED0 half-open peak, same predictions",
               "repeats": "first three admitted rows, same process after original pass; no outcome selection",
               "purpose": "paired development prediction table; no accuracy claim" if full_development else
                          "geometry/numerical/cost feasibility only; no accuracy claim",
               "invariant_checks": "passed"}
    checkpoint_root = C.CHECKPOINT.parent
    receipt["checkpoint_receipts"] = {name: C.sha256_file(checkpoint_root/name)
                                      for name in ["acquisition.json", "checkpoint_inventory.tsv"]}
    save(out / "receipt.json", receipt)
    (out / "executed_source.py").write_bytes(Path(__file__).read_bytes())
    fasta = pysam.FastaFile(C.FASTA_PATH)
    accepted, rejected = [], []
    for row in rows.to_dict("records"):
        ref, alt = row["ref"], row["alt"]
        if len(ref) != 1 or len(alt) != 1 or ref not in C.ACGT or alt not in C.ACGT or ref == alt:
            raise C.ContractError("Manifest has non-SNV or equal alleles.")
        lo, hi = C.window_bounds(int(row["pos_hg38"]), 1048576)
        if (lo, hi) != (int(row["window_start0"]), int(row["window_end0"])):
            raise C.ContractError("Manifest window differs from native window.")
        check = C.validate_window(fasta, row["chr"], lo, hi, int(row["pos_hg38"]), ref)
        if not check["acgt_ok"] or not check["reference_match"]:
            rejected.append({"key": row["key"], "geometry_stratum": row["geometry_stratum"], **check})
        else:
            accepted.append(row)
    save(out / "sequence_admission.json", {"admitted": len(accepted), "rejected": rejected,
                                           "replacements": 0, "outcomes_read": False})
    if not accepted:
        raise C.ContractError("No sequence-admitted fixed cases.")
    restore_tick = time.monotonic()
    model, device, load_s = C.load_model()
    jax.block_until_ready((model._params, model._state))
    receipt["actual_device"] = {"kind": str(device.device_kind), "id": int(device.id),
                                "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
                                "synchronized_restore_seconds": time.monotonic()-restore_tick}
    if "L40S" not in str(device.device_kind).upper():
        raise C.ContractError("JAX's actual allocated device is not L40S.")
    save(out / "receipt.json", receipt)
    outputs = [dna_client.OutputType.ATAC, dna_client.OutputType.DNASE]

    def predict(row):
        tick = time.monotonic()
        lo, hi = int(row["window_start0"]), int(row["window_end0"])
        seqs = C.extract_ref_alt(fasta, row["chr"], lo, hi, int(row["pos_hg38"]), row["ref"], row["alt"])
        interval = genome.Interval(chromosome=row["chr"], start=lo, end=hi)
        variant = genome.Variant(chromosome=row["chr"], position=int(row["pos_hg38"]),
                                 reference_bases=row["ref"], alternate_bases=row["alt"])
        tracks = [model.predict_sequence(sequence=seq, requested_outputs=outputs,
                  ontology_terms=C.LIVER_TERMS, interval=interval).atac for seq in seqs]
        names = None
        values = []
        for td in tracks:
            if td is None or td.interval is None or td.resolution != 1 or td.interval.chromosome != interval.chromosome or \
                    (td.interval.start, td.interval.end) != (lo, hi):
                raise C.ContractError("Returned ATAC geometry differs; no clipping/cropping accepted.")
            v = np.asarray(td.values, dtype=np.float32)
            current_names = td.metadata["name"].astype(str).to_numpy()
            if v.ndim != 2 or v.shape != (hi-lo, len(current_names)) or not np.isfinite(v).all() or (v < 0).any():
                raise C.ContractError("Invalid returned ATAC shape or scale.")
            if names is not None and not np.array_equal(names, current_names):
                raise C.ContractError("REF/ALT track identities differ.")
            names = current_names
            values.append(jnp.asarray(v))
        indices = []
        for name in C.AG_ATAC_TRACKS:
            hits = np.flatnonzero(names == name)
            if len(hits) != 1:
                raise C.ContractError("Expected exactly one named track: " + name)
            indices.append(int(hits[0]))
        masks = {"center": make_mask(interval.as_unstranded(), variant, width=501, resolution=1),
                 "peak": peak_mask(tracks[0].interval, int(row["peak_start0"]), int(row["peak_end0"]))}
        center_indices = np.flatnonzero(masks["center"][:, 0])
        if len(center_indices) != 501 or (lo+int(center_indices[0]), lo+int(center_indices[-1])+1) != \
                (int(row["mask_start0"]), int(row["mask_end0"])):
            raise C.ContractError("Installed center mask differs from geometry census.")
        result = {"key": row["key"], "peak_id": row["peak_id"], "heldout_fold": row["heldout_fold"],
                  "geometry_stratum": row["geometry_stratum"], "peak_width_bp": int(masks["peak"].sum()),
                  "sequence_sha256": hashlib.sha256(seqs[0].encode()).hexdigest()}
        for label, mask in masks.items():
            score = np.asarray(aggregate(*values, jnp.asarray(mask), aggregation_type=scorer.aggregation_type),
                               dtype=np.float64)[indices]
            if not np.isfinite(score).all():
                raise C.ContractError("Nonfinite native readout.")
            result[label + "_score"] = float(score.mean())
            for j, value in enumerate(score):
                result[f"{label}_track{j}"] = float(value)
        result["seconds"] = time.monotonic()-tick
        return result

    results = []
    with (out / "scores.tsv").open("x", newline="") as handle:
        writer = None
        for row in accepted:
            result = predict(row)
            if writer is None:
                writer = csv.DictWriter(handle, fieldnames=list(result), delimiter="\t")
                writer.writeheader()
            writer.writerow(result)
            handle.flush()
            results.append(result)
            print(json.dumps({"scored": len(results), "total_admitted": len(accepted),
                              "seconds": result["seconds"]}), flush=True)
    repeats = []
    for row, original in zip(accepted[:3], results[:3]):
        repeat = predict(row)
        repeats.append({"key": row["key"], "center_delta": repeat["center_score"]-original["center_score"],
                        "peak_delta": repeat["peak_score"]-original["peak_score"], "repeat": repeat})
    warm = [row["seconds"] for row in results[1:]]
    save(out / "summary.json", {"status": "development_paired_readouts_complete" if full_development else
                                "geometry_feasibility_complete", "selected": len(rows),
         "admitted": len(accepted), "scored": len(results), "rejected": len(rejected),
         "load_seconds": load_s, "first_pair_seconds": results[0]["seconds"],
         "median_warm_pair_seconds": float(np.median(warm)) if warm else None,
         "same_process_repeats": repeats, "device_memory": C.memory_stats(device),
         "outcomes_read": False, "accuracy_evaluated": False})
    fasta.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--full-development", action="store_true",
                        help="Use the fixed 25,477-row development manifest; fold 0 remains excluded.")
    args = parser.parse_args()
    main(args.manifest, args.full_development)
