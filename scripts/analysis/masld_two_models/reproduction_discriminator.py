#!/usr/bin/env python3
"""Distinguish changed inference geometry from cross-process numerical drift.

The unmodified native scoring script runs twice in separate processes. Each
process scores the same 32 metadata-only variants twice. Saved adapter weights
are evaluated at batch one and at their original batch-four companion/order
geometry, twice per process and in two processes per geometry. No association
outcomes, fitting, precision changes, or revised tolerance are involved.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import pandas as pd

import variant_inference as V

FOLD = 1
N = 32
ATOL = RTOL = 1e-4
PRODUCER = V.ROOT/"scripts/analysis/alphagenome_program/i1_01_score_local.py"
NATIVE = V.BASE/"native/native_1048576.tsv"
SHORT = V.BASE/"native_lengths_21775189/native_2048.tsv"
INVENTORY = V.ROOT/"Analysis/MASLD_Model_Benchmark/executions/alphagenome-weights-20260915T103442Z/checkpoint_inventory.tsv"


def sha256(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4*1024*1024), b""):
            value.update(block)
    return value.hexdigest()


def checkpoint_receipt():
    rows = list(csv.DictReader(INVENTORY.open(), delimiter="\t"))
    for row in rows:
        path = V.C.CHECKPOINT/row["path"]
        if path.stat().st_size != int(row["bytes"]) or sha256(path) != row["sha256"]:
            raise ValueError("Checkpoint differs from acquired inventory: " + row["path"])
    return {"inventory_sha256": sha256(INVENTORY), "files_verified": len(rows),
            "bytes_verified": sum(int(row["bytes"]) for row in rows),
            "checkpoint_unchanged_from_acquisition": True}


def source_adapters():
    folder = V.BASE/f"nested_f{FOLD}_{21783074+FOLD}"
    records = {}
    ordering = None
    for seed in V.SEEDS:
        frame = pd.read_csv(folder/f"{V.RECIPE}__seed{seed}"/"validation_predictions.tsv", sep="\t",
                            usecols=["variant_id", "validation_fold", "predicted_beta"])
        frame["key"] = frame.variant_id.str.removeprefix("chr")
        if frame.key.duplicated().any() or not frame.validation_fold.eq(FOLD).all():
            raise ValueError("Saved adapter identity or fold differs")
        keys = frame.key.to_numpy()
        if ordering is not None and not np.array_equal(ordering, keys):
            raise ValueError("Saved seed prediction order differs")
        ordering = keys
        records[seed] = frame.set_index("key").predicted_beta
    return ordering, records


def load_manifest(path):
    frame = pd.read_csv(path, sep="\t")
    expected = {"key", "chr", "pos_hg38", "ref", "alt", "heldout_fold", "selection_hash"}
    if set(frame) != expected or len(frame) != N or frame.key.duplicated().any() or not frame.heldout_fold.eq(FOLD).all():
        raise ValueError("Expected the same 32 outcome-blind fold-1 variants")
    for row in frame.to_dict("records"):
        if row["key"] != V.C.key_of(row["chr"], row["pos_hg38"], row["ref"], row["alt"]):
            raise ValueError("Manifest identity disagrees with coordinates")
    return frame


def native_manifest(frame):
    metadata = pd.read_csv(NATIVE, sep="\t", usecols=["key", "block_1mb"]).set_index("key")
    result = frame.drop(columns="selection_hash").copy()
    result["block_1mb"] = metadata.loc[result.key, "block_1mb"].to_numpy()
    # i1_01_score_local copies this field to its output and never uses it for
    # inference. Every value is synthetic zero, never a measured molecular beta.
    result["beta_alt"] = 0.0
    result = pd.concat([result.assign(repeat_index=repeat) for repeat in range(2)], ignore_index=True)
    return result


def adapter_child(args):
    import jax
    import jax.numpy as jnp
    import model_scalar as M
    if args.out.exists():
        raise FileExistsError(args.out)
    frame = load_manifest(args.manifest)
    order, saved = source_adapters()
    metadata = pd.read_csv(SHORT, sep="\t", usecols=["key", "chr", "pos_hg38", "ref", "alt", "heldout_fold"])
    metadata = metadata.set_index("key")
    if not set(order).issubset(metadata.index) or not metadata.loc[order].heldout_fold.eq(FOLD).all():
        raise ValueError("Original adapter batch companions lack eligible metadata")
    position = {key: index for index, key in enumerate(order)}
    predictor = V.AdapterPredictor(FOLD)
    paths = M.qv_paths(predictor.base, 5)

    # Match model_train.py's inference function including its zero native-score
    # input, frozen state and zero native slope. No optimizer or label is used.
    def original_predict(weights, backbone, stats, ref, alt, pool, native):
        return M.sequence_effect(predictor.trunk, backbone, stats, weights, paths,
                                 "adapter", ref, alt, pool, native, 0.0)

    predict = jax.jit(original_predict)
    rows = []
    for repeat in range(2):
        for selected in frame.to_dict("records"):
            key = selected["key"]
            if args.batch_size == 1:
                companions = [key]
                selected_offset = 0
            else:
                first = (position[key]//4)*4
                companions = list(order[first:first+4])
                companions += [companions[-1]]*(4-len(companions))
                selected_offset = companions.index(key)
            sequences = [V.window(predictor.fasta, metadata.loc[companion].to_dict(), 2048)[0]
                         for companion in companions]
            ref = jnp.asarray(np.stack([predictor.encoder.encode(sequence[0]) for sequence in sequences]), dtype=jnp.float32)
            alt = jnp.asarray(np.stack([predictor.encoder.encode(sequence[1]) for sequence in sequences]), dtype=jnp.float32)
            pool = jnp.repeat(predictor.pool, args.batch_size, axis=0)
            native = jnp.zeros(args.batch_size, dtype=jnp.float32)
            for seed, weights in zip(V.SEEDS, predictor.weights):
                tick = time.monotonic()
                values = np.asarray(predict(weights, predictor.base, predictor.state, ref, alt, pool, native))
                value = float(values[selected_offset])
                if not np.isfinite(value):
                    raise ValueError("Nonfinite adapter prediction")
                rows.append({"key": key, "seed": seed, "repeat_index": repeat,
                             "process_repeat": args.process_repeat, "batch_size": args.batch_size,
                             "batch_companions": ";".join(companions), "selected_batch_offset": selected_offset,
                             "prediction": value, "saved_prediction": float(saved[seed].loc[key]),
                             "prediction_minus_saved": value-float(saved[seed].loc[key]),
                             "synchronized_call_seconds": time.monotonic()-tick})
        # Preserve completed repeats if a subsequent numerical or runtime issue occurs.
        pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    from inference_cost import gpu_identity
    receipt = {"batch_size": args.batch_size, "process_repeat": args.process_repeat,
               "variant_queries": N, "within_process_repeats": 2, "seeds": list(V.SEEDS),
               "card": predictor.card, "gpu_processes": gpu_identity(), "host": platform.node(),
               "jax": jax.__version__, "numpy": np.__version__,
               "weights_sha256": predictor.weight_sha256, "memory": V.C.memory_stats(predictor.device),
               "measured_outcomes_loaded": False, "precision_changed": False,
               "batch_companions": "same order and padding as saved validation batch when batch_size=4"}
    args.out.with_suffix(".json").write_text(json.dumps(receipt, indent=2)+"\n")


def comparison(a, b, label):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("Invalid paired comparison: " + label)
    difference = a-b
    return {"comparison": label, "n": a.size, "max_absolute_difference": float(np.max(np.abs(difference))),
            "RMS_difference": float(np.sqrt(np.mean(difference*difference))),
            "exactly_identical": bool(np.array_equal(a, b)),
            "passes_original_1e_minus4_tolerance": bool(np.allclose(a, b, rtol=RTOL, atol=ATOL)),
            "prediction_sign_changes": int(np.sum(np.sign(a) != np.sign(b))),
            "atol": ATOL, "rtol": RTOL}


def summarize_native(directory, frame):
    source = pd.read_csv(NATIVE, sep="\t", usecols=["key", "local_atac_liver"]).set_index("key")
    comparisons, rows, first = [], [], []
    expected_keys = np.tile(frame.key.to_numpy(), 2)
    for process in range(2):
        observed = pd.read_csv(directory/f"native_process{process}.tsv", sep="\t")
        if not np.array_equal(observed.key.to_numpy(), expected_keys) or not observed.beta_alt.eq(0).all():
            raise ValueError("Original producer changed query order or synthetic placeholder")
        observed["repeat_index"] = observed.groupby("key", sort=False).cumcount()
        observed["process_repeat"] = process
        observed["saved_prediction"] = source.loc[observed.key, "local_atac_liver"].to_numpy()
        observed["prediction_minus_saved"] = observed.local_atac_liver-observed.saved_prediction
        rows.append(observed.drop(columns="beta_alt"))
        a = observed.loc[observed.repeat_index.eq(0), "local_atac_liver"].to_numpy()
        b = observed.loc[observed.repeat_index.eq(1), "local_atac_liver"].to_numpy()
        first.append(a)
        comparisons.append(comparison(a, b, f"native_within_process_{process}"))
        comparisons.append(comparison(a, source.loc[frame.key, "local_atac_liver"], f"native_process_{process}_vs_saved"))
    comparisons.append(comparison(first[0], first[1], "native_between_processes"))
    pd.concat(rows, ignore_index=True).to_csv(directory/"native_per_site_diagnostics.tsv", sep="\t", index=False)
    return comparisons


def summarize_adapter(directory):
    tables, summaries, ensembles = {}, [], {}
    for batch in (1, 4):
        for process in range(2):
            frame = pd.read_csv(directory/f"adapter_batch{batch}_process{process}.tsv", sep="\t")
            if len(frame) != N*5*2:
                raise ValueError("Adapter diagnostic census differs")
            a = frame.loc[frame.repeat_index.eq(0)].set_index(["key", "seed"]).sort_index()
            b = frame.loc[frame.repeat_index.eq(1)].set_index(["key", "seed"]).sort_index()
            if not a.index.equals(b.index):
                raise ValueError("Adapter within-process identity differs")
            tables[batch, process] = a
            ensemble = a.groupby("key")[["prediction", "saved_prediction"]].mean().sort_index()
            ensembles[batch, process] = ensemble
            summaries.append(comparison(a.prediction, b.prediction, f"adapter_batch{batch}_within_process{process}_seedwise"))
            summaries.append(comparison(a.prediction, a.saved_prediction, f"adapter_batch{batch}_process{process}_vs_saved_seedwise"))
            summaries.append(comparison(ensemble.prediction, ensemble.saved_prediction,
                                        f"adapter_batch{batch}_process{process}_vs_saved_ensemble"))
        summaries.append(comparison(tables[batch, 0].prediction, tables[batch, 1].prediction,
                                    f"adapter_batch{batch}_between_processes_seedwise"))
        summaries.append(comparison(ensembles[batch, 0].prediction, ensembles[batch, 1].prediction,
                                    f"adapter_batch{batch}_between_processes_ensemble"))
    for process in range(2):
        summaries.append(comparison(ensembles[1, process].prediction, ensembles[4, process].prediction,
                                    f"adapter_batch1_vs_batch4_process_pair{process}_ensemble"))
    return summaries


def main(args):
    V.require_gpu()
    if args.arm == "adapter":
        adapter_child(args)
        return
    if args.out.exists():
        raise FileExistsError(args.out)
    args.out.mkdir(parents=True)
    frame = load_manifest(args.manifest)
    checkpoint = checkpoint_receipt()
    design = {"measured_outcomes_loaded": False, "source_manifest": str(args.manifest),
              "source_manifest_sha256": sha256(args.manifest), "fold": FOLD, "variants": N,
              "native_placeholder": "beta_alt=0 is an unused synthetic placeholder required by the unmodified scorer's output schema; it is not a measured effect",
              "native": "unchanged original producer, two independent processes, same32 queries repeated twice per process",
              "adapter": "saved weights, original model_train inference function, batch1 and original validation batch4 companions, two processes per geometry and two query repeats per process",
              "tolerance_unchanged": {"atol": ATOL, "rtol": RTOL}, "precision_changed": False,
              "limits": "diagnostic predictions only; no model accuracy or speed claim; batch geometry and process both vary between adapter geometry arms"}
    (args.out/"design.json").write_text(json.dumps(design, indent=2)+"\n")
    source = native_manifest(frame)
    native_input = args.out/"native_metadata_with_synthetic_unused_placeholder.tsv"
    source.to_csv(native_input, sep="\t", index=False)
    code_paths = [Path(__file__), Path(V.__file__), PRODUCER,
                  V.ROOT/"scripts/analysis/alphagenome_program/i1_common.py",
                  V.ROOT/"scripts/analysis/alphagenome_campaign/model_scalar.py",
                  V.ROOT/"scripts/analysis/alphagenome_campaign/model_train.py",
                  V.ROOT/"scripts/analysis/alphagenome_program/i1_extract_mpra_embeddings.py"]
    code_hashes = {str(path): sha256(path) for path in code_paths}
    for process in range(2):
        out = args.out/f"native_process{process}.tsv"
        with (args.out/f"native_process{process}.log").open("w") as log:
            subprocess.run([sys.executable, str(PRODUCER), "--variants", str(native_input), "--length", "1048576",
                            "--out", str(out), "--require-card", "l40s", "--flush-every", "1"],
                           check=True, stdout=log, stderr=subprocess.STDOUT)
    native_summary = summarize_native(args.out, frame)
    (args.out/"native_comparisons.json").write_text(json.dumps(native_summary, indent=2)+"\n")
    for batch in (1, 4):
        for process in range(2):
            out = args.out/f"adapter_batch{batch}_process{process}.tsv"
            with (args.out/f"adapter_batch{batch}_process{process}.log").open("w") as log:
                subprocess.run([sys.executable, str(Path(__file__).resolve()), "--arm", "adapter",
                                "--manifest", str(args.manifest), "--batch-size", str(batch),
                                "--process-repeat", str(process), "--out", str(out)],
                               check=True, stdout=log, stderr=subprocess.STDOUT)
    adapter_summary = summarize_adapter(args.out)
    for path, digest in code_hashes.items():
        if sha256(path) != digest:
            raise ValueError("Code changed during diagnostic: " + path)
    result = {"state": "prediction-only reproduction discriminator; fixed weights and unchanged tolerance",
              "native": native_summary, "adapter": adapter_summary, "checkpoint": checkpoint,
              "code_sha256": code_hashes, "manifest_sha256": sha256(args.manifest),
              "measured_outcomes_loaded": False, "fit_performed": False,
              "environment": {"host": platform.node(), "python": sys.version,
                              "numpy": np.__version__, "slurm_job_id": os.environ["SLURM_JOB_ID"]}}
    (args.out/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--arm", choices=["adapter"])
    parser.add_argument("--batch-size", type=int, choices=[1, 4])
    parser.add_argument("--process-repeat", type=int, choices=[0, 1])
    main(parser.parse_args())
