#!/usr/bin/env python3
"""Same-node, same-GPU, batch-one timing without molecular outcome fields.

Run the native and five-seed adapter arms in separate child processes so GPU
peak memory and cold initialization are attributable to one arm. No fitting,
calibration, hyperparameter selection, or measured-effect evaluation occurs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time

import numpy as np
import pandas as pd

import variant_inference as V

NATIVE = V.BASE/"native/native_1048576.tsv"
FOLD = 1
N_VARIANTS = 32
WARMUPS = 3


def manifest():
    columns = ["key", "chr", "pos_hg38", "ref", "alt", "heldout_fold"]
    native = pd.read_csv(NATIVE, sep="\t", usecols=columns)
    native = native.loc[native.heldout_fold.eq(FOLD)].copy()
    if native.key.duplicated().any():
        raise ValueError("Duplicate native variant identity")
    shared = set(native.key)
    seed_keys = []
    folder = V.BASE/f"nested_f{FOLD}_{21783074+FOLD}"
    for seed in V.SEEDS:
        frame = pd.read_csv(folder/f"{V.RECIPE}__seed{seed}"/"validation_predictions.tsv",
                            sep="\t", usecols=["variant_id", "validation_fold"])
        if frame.variant_id.duplicated().any() or not frame.validation_fold.eq(FOLD).all():
            raise ValueError("Adapter identity or development fold differs")
        seed_keys.append(set(frame.variant_id.str.removeprefix("chr")))
        shared &= seed_keys[-1]
    if any(keys != seed_keys[0] for keys in seed_keys[1:]):
        raise ValueError("Adapter seed coverage differs")
    frame = native.loc[native.key.isin(shared)].copy()
    frame["selection_hash"] = frame.key.map(lambda value: hashlib.sha256(
        ("masld-inference-cost-20260923|"+value).encode()).hexdigest())
    frame = frame.sort_values(["selection_hash", "key"]).head(N_VARIANTS)
    if len(frame) != N_VARIANTS:
        raise ValueError("Too few shared development variants")
    return frame.reset_index(drop=True)


def gpu_identity():
    # Query the GPU actually hosting this process, rather than the first GPU
    # listed on a multi-GPU node outside CUDA_VISIBLE_DEVICES.
    result = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,gpu_uuid,used_memory",
                             "--format=csv,noheader,nounits"], check=True, capture_output=True, text=True)
    rows = []
    for line in result.stdout.splitlines():
        parts = [part.strip() for part in line.split(",")]
        if len(parts) == 3 and parts[0] == str(os.getpid()):
            rows.append({"gpu_uuid": parts[1], "observed_process_memory_mib": parts[2]})
    if not rows:
        raise RuntimeError("Could not identify the GPU hosting the timed process")
    return rows


def verify_saved_predictions(arm, frame, values):
    # Only predictions and identity are read. Beta/association outcomes remain
    # excluded even in the numerical reproduction check.
    if arm.startswith("native1m"):
        source = pd.read_csv(NATIVE, sep="\t", usecols=["key", "local_atac_liver"]).set_index("key")
        expected = source.loc[frame.key, "local_atac_liver"].to_numpy(float)
    else:
        folder = V.BASE/f"nested_f{FOLD}_{21783074+FOLD}"
        arrays = []
        for seed in V.SEEDS:
            source = pd.read_csv(folder/f"{V.RECIPE}__seed{seed}"/"validation_predictions.tsv",
                                 sep="\t", usecols=["variant_id", "predicted_beta"])
            source["key"] = source.variant_id.str.removeprefix("chr")
            arrays.append(source.set_index("key").loc[frame.key, "predicted_beta"].to_numpy(float))
        expected = np.mean(arrays, axis=0)
    maximum = float(np.max(np.abs(values-expected)))
    return {"passed": bool(np.allclose(values, expected, rtol=1e-4, atol=1e-4)),
            "max_absolute_difference": maximum, "atol": 1e-4, "rtol": 1e-4,
            "comparison": "predictions only; no measured effect used",
            "per_variant": [{"key": key, "saved_prediction": float(saved),
                             "absolute_difference": float(abs(value-saved))}
                            for key, value, saved in zip(frame.key, values, expected)]}


def child(args):
    if args.out.exists():
        raise FileExistsError(args.out)
    frame = pd.read_csv(args.manifest, sep="\t")
    forbidden = {"beta_alt", "observed_beta", "q_val", "p_nominal", "p_val"}
    if len(frame) != N_VARIANTS or forbidden.intersection(frame.columns) or not frame.heldout_fold.eq(FOLD).all():
        raise ValueError("Timing manifest population or columns differ")
    genome_stat = Path(V.C.FASTA_PATH).stat()
    genome_identity = {"path": V.C.FASTA_PATH, "size": genome_stat.st_size,
                       "mtime_ns": genome_stat.st_mtime_ns,
                       "fai_sha256": hashlib.sha256(Path(V.C.FASTA_PATH+".fai").read_bytes()).hexdigest()}
    predictor = (V.NativePredictor(atac_only=args.arm == "native1m_atac_only")
                 if args.arm.startswith("native1m") else V.AdapterPredictor(FOLD))
    warmup = []
    for index in range(WARMUPS):
        tick = time.monotonic()
        predictor.predict(frame.iloc[index].to_dict())
        warmup.append(time.monotonic()-tick)
    timings, values = [], []
    for _, row in frame.iterrows():
        tick = time.monotonic()
        values.append(predictor.predict(row.to_dict()))
        timings.append(time.monotonic()-tick)
    # Persist actual calls before any comparison can fail, so discrepancies
    # remain inspectable without rerunning a molecular-effect evaluation.
    calls = [{"key": key, "seconds": sec, "prediction": value}
             for key, sec, value in zip(frame.key, timings, values)]
    pd.DataFrame(calls).to_csv(args.out.with_suffix(".calls.tsv"), sep="\t", index=False)
    verification = verify_saved_predictions(args.arm, frame, np.asarray(values))
    import jax
    memory = V.C.memory_stats(predictor.device)
    stat_after = Path(V.C.FASTA_PATH).stat()
    if (stat_after.st_size, stat_after.st_mtime_ns) != (genome_stat.st_size, genome_stat.st_mtime_ns):
        raise ValueError("Reference genome changed during inference")
    weight_hashes = getattr(predictor, "weight_sha256", {})
    for path, digest in weight_hashes.items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError("Adapter weights changed during inference")
    result = {"arm": args.arm, "variants": N_VARIANTS, "variant_batch_size": 1,
              "alleles_per_variant": 2, "adapter_seeds": 5 if args.arm == "adapter" else None,
              "warmup_calls": WARMUPS, "initialization_seconds": predictor.initialization_seconds,
              "checkpoint_restore_seconds": predictor.restore_seconds,
              "first_call_seconds_including_compile": warmup[0],
              "later_warmup_seconds": warmup[1:], "warm_seconds_mean": float(np.mean(timings)),
              "warm_seconds_median": float(np.median(timings)),
              "warm_seconds_q10_q90": np.quantile(timings, [.1, .9]).tolist(),
              "timing_scope": "FASTA validation/extraction, encoding, both alleles, readout, host-synchronized scalar; adapter includes all five seeds",
              "native_requested_outputs": getattr(predictor, "requested_output_names", None),
              "jax_memory_stats": memory, "host_peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              "gpu_processes": gpu_identity(), "card": predictor.card,
              "host": platform.node(), "python": sys.version, "jax": jax.__version__, "numpy": np.__version__,
              "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
              "genome": genome_identity, "checkpoint": str(V.C.CHECKPOINT),
              "adapter_weights_sha256": weight_hashes, "prediction_reproduction": verification,
              "measured_outcomes_loaded": False, "weights_changed": False,
              "slurm_job_id": os.environ["SLURM_JOB_ID"],
              "per_variant": calls}
    args.out.write_text(json.dumps(result, indent=2)+"\n")
    if not verification["passed"] and args.arm != "native1m_atac_only":
        raise ValueError(f"{args.arm} differs from saved predictions; retained values and timing in {args.out}")


def main(args):
    V.require_gpu()
    if args.arm:
        child(args)
        return
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.mkdir(parents=True)
    frame = manifest()
    path = args.out/"outcome_blind_variants.tsv"
    frame.to_csv(path, sep="\t", index=False)
    arms = {}
    for arm in ("native1m", "native1m_atac_only", "adapter"):
        destination = args.out/(arm+".json")
        with (args.out/(arm+".log")).open("w") as log:
            call = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--arm", arm,
                                   "--manifest", str(path), "--out", str(destination)],
                                  check=False, stdout=log, stderr=subprocess.STDOUT)
        if not destination.exists():
            raise RuntimeError(f"{arm} failed before completing timing (exit {call.returncode}); inspect {arm}.log")
        arms[arm] = json.loads(destination.read_text())
    expected_gpu = {item["gpu_uuid"] for item in arms["native1m"]["gpu_processes"]}
    if any(record["host"] != arms["native1m"]["host"] or
           {item["gpu_uuid"] for item in record["gpu_processes"]} != expected_gpu
           for record in arms.values()):
        raise ValueError("Native and adapter timing used different hardware")
    valid = all(arms[arm]["prediction_reproduction"]["passed"] for arm in ("native1m", "adapter"))
    result = {"state": "candidate inference-only cost; not a noninferiority analysis",
              "variants": N_VARIANTS, "development_fold": FOLD, "selection": "SHA256 of identity only",
              "measured_outcomes_loaded": False, "same_hardware_verified": True,
              "arm_order": ["native1m", "native1m_atac_only", "adapter"], "order_counterbalanced": False,
              "primary_prediction_reproduction_passed": valid,
              "variant_batch_size": 1, "warmup_calls_per_arm": WARMUPS,
              "warm_seconds_mean": {arm: arms[arm]["warm_seconds_mean"] for arm in ("native1m", "adapter")},
              "native_to_five_seed_adapter_mean_time_ratio": arms["native1m"]["warm_seconds_mean"]/arms["adapter"]["warm_seconds_mean"] if valid else None,
              "combined_warm_seconds_estimated_sum": sum(arms[arm]["warm_seconds_mean"] for arm in ("native1m", "adapter")) if valid else None,
              "combined_cost_label": "estimate from separately measured model calls; excludes linear combination and process startup",
              "native_primary_outputs": "ATAC+DNase, matching the saved evaluated scoring graph; only ATAC is read out",
              "ATAC_only_diagnostic": {"reproduces_saved_predictions": arms["native1m_atac_only"]["prediction_reproduction"]["passed"],
                                       "max_absolute_difference_from_saved": arms["native1m_atac_only"]["prediction_reproduction"]["max_absolute_difference"],
                                       "warm_seconds_mean": arms["native1m_atac_only"]["warm_seconds_mean"],
                                       "role": "output-request diagnostic; not used for the primary cost ratio"},
              "limits": "32 selected eligible variants, one L40S, batch-one latency; optimized throughput and hardware-to-hardware variation not measured",
              "manifest_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "script_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (Path(__file__), Path(V.__file__))},
              "arms": {arm: arm+".json" for arm in arms}}
    (args.out/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)
    if not valid:
        raise RuntimeError("Primary prediction reproduction failed; retained diagnostics, no admitted cost ratio")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--arm", choices=["native1m", "adapter", "native1m_atac_only"])
    parser.add_argument("--manifest", type=Path)
    main(parser.parse_args())
