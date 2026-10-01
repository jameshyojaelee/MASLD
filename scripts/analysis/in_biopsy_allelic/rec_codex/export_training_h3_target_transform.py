"""Reconstruct only the frozen all-99 H3 concentration target transform.

Protocol v1: no RNA values, receiving data, prediction fitting, or target
selection. Use exact source functions extracted without importing their
top-level data loaders. Every failure stops before coefficient export.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
BENCH = ROOT / "Analysis/MASLD_Model_Benchmark"
FIX = BENCH / "executions/model-data-064-21079902/fixture/molecular"
STABLE = BENCH / "executions/chromatin-stable-rrr-forms-20260908T192024Z"
RELEASE = BENCH / "release/masld-liver-chromatin-state-v1.2"
SOURCE = {
    STABLE / "transfer_functions.py": "05b995e349fb0b14a13201310523dd7428ef8a4009190622fba7a31c5cb7e117",
    STABLE / "run_stable.py": "0ebec6e0a319879af8fec482ccf06d2750cb47e0341d3015102a93350fb1c25b",
    BENCH / "executions/chromatin-beyond-histology-20260907T161521Z/common.py": "f2a285bcbf0197a252bf01faa268a113b7547fca3705125e46005bb63ea7ce3e",
    BENCH / "executions/chromatin-state-v12-build-20260908T172230Z/build_v12.py": "aacedf36446e5f5bdeacbce08cf03cf55030fa3dcc1c79d6ad7226b25395c77d",
}
INPUTS = {
    FIX / "h3k27ac_counts.npy": "9a4fa736a6dd5de0cd44e51428f2bf19e8dae9adb7a279e5b6c8f1c14103cd22",
    FIX / "h3k27ac_feature_axis.tsv": "cbe35aeb188e531283a1bcd636ba14de4cfbf215dfc544b99e3ebb97f9831986",
    FIX / "participant_axis.tsv": "7e5be1a16df505e5daad451398e8b1a4b3cbbc48a32cd1381113e1ae4b7347bc",
    FIX / "receipt.json": "5cda57f2de6f9586d392ee3a2c82dc74e00835a1029006a0df844009a6844317",
    RELEASE / "weights/chromatin_state_v1_2.npz": "72f9eaba359f7c4c61cac2e00da2d5539b5d3ec868dc33b1bacc2ef0105694ec",
    STABLE / "out/stable_heads_all99.npz": "fa3e8b6e82d672e1df5c10d0abb1fe033f4042381c8ed62135647069df8ebe28",
}
DESCRIPTORS = ("intercept", "gini_region_counts", "shannon_entropy_natural_log",
               "iqr_log2_positive_counts", "top5pct_region_count_fraction")
# All sources and model means are float64; 1e-9 log2-CPM accommodates BLAS
# least-squares roundoff, without admitting a float32 target approximation.
ABS_TOL = 1e-9


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_axis(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def source_function(path, name, namespace):
    text = path.read_text()
    nodes = [n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name == name]
    require(len(nodes) == 1, "Ambiguous frozen source function " + name)
    source = ast.get_source_segment(text, nodes[0])
    exec(compile(source, str(path) + "::" + name, "exec"), namespace)
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Read development counts only inside a compute job")
    require(not args.output.exists(), "Output exists; refusing overwrite")
    hashes = {}
    for path, expected in {**SOURCE, **INPUTS}.items():
        actual = sha(path)
        require(actual == expected, "Frozen source/input hash changed: " + str(path))
        hashes[str(path)] = actual
    participant_rows = read_axis(FIX / "participant_axis.tsv")
    region_rows = read_axis(FIX / "h3k27ac_feature_axis.tsv")
    ids = np.array([r["participant_id"] for r in participant_rows], dtype=str)
    keys = np.array([r["opaque_source_feature_key"] for r in region_rows], dtype=str)
    require(len(ids) == len(set(ids)) == 99, "Training participant roster changed")
    require(len(keys) == len(set(keys)) == 96460, "Complete region axis changed")
    require([int(r["participant_index"]) for r in participant_rows] == list(range(99)), "Participant order changed")
    require([int(r["h3k27ac_feature_index"]) for r in region_rows] == list(range(96460)), "Region order changed")
    receipt = json.loads((FIX / "receipt.json").read_text())
    require(receipt["h3k27ac"]["shape"] == [99, 96460] and receipt["participants"] == 99,
            "Fixture receipt does not identify complete training target")
    with np.load(RELEASE / "weights/chromatin_state_v1_2.npz", allow_pickle=True) as weights:
        require(int(weights["n_train"]) == int(weights["prof_n_train"]) == 99, "Released training n differs")
        require(np.array_equal(weights["prof_region_index"], np.arange(96460)), "Released complete region index differs")
        require(np.array_equal(weights["prof_region_key"].astype(str), keys), "Released region identity differs")
        require(len(weights["gene_axis"]) == 42163, "Released RNA axis size differs")
        mean_release = np.asarray(weights["prof_rrr_ym"], dtype=np.float64).copy()
    with np.load(STABLE / "out/stable_heads_all99.npz", allow_pickle=False) as heads:
        mean_stable = np.asarray(heads["rrr_ym"], dtype=np.float64).copy()
    require(np.array_equal(mean_release, mean_stable), "Released and final all-99 target means differ")
    require(mean_release.shape == (96460,) and np.isfinite(mean_release).all(), "Invalid stored target mean")
    namespace = {"np": np}
    sources = {}
    for name in ("conc_block", "residualise"):
        sources[name] = source_function(STABLE / "transfer_functions.py", name, namespace)
    sources["logcpm"] = source_function(BENCH / "executions/chromatin-beyond-histology-20260907T161521Z/common.py", "logcpm", namespace)
    # This is the sole molecular value read: allowed all-99 development H3.
    native = np.load(FIX / "h3k27ac_counts.npy", mmap_mode="r", allow_pickle=False)
    require(native.shape == (99, 96460) and native.dtype == np.uint32, "Count shape/dtype changed")
    counts = np.asarray(native, dtype=np.float64)
    require(np.isfinite(counts).all() and (counts >= 0).all() and (counts.sum(1) > 0).all(), "Invalid development H3 counts")
    require((counts > 0).any(1).all(), "No positive counts for a training donor")
    logcpm = namespace["logcpm"](counts)
    descriptors = namespace["conc_block"](counts)
    require(descriptors.shape == (99, 4) and np.isfinite(descriptors).all(), "Invalid source descriptors")
    design = np.column_stack([np.ones(99), descriptors])
    require(np.linalg.matrix_rank(design) == 5, "Concentration design is not rank five")
    coef, _, rank, singular = np.linalg.lstsq(design, logcpm, rcond=None)
    require(rank == 5 and coef.shape == (5, 96460) and np.isfinite(coef).all(), "Invalid target coefficients")
    reconstructed = logcpm - design @ coef
    original = namespace["residualise"](logcpm, descriptors, np.arange(99))
    source_difference = float(np.max(np.abs(reconstructed - original)))
    mean_difference = float(np.max(np.abs(reconstructed.mean(0) - mean_release)))
    require(source_difference <= ABS_TOL, "Original source all-99 residual convention not reproduced")
    require(mean_difference <= ABS_TOL, "Reconstructed target mean differs from released all-99 mean")
    # Reversible library scaling must preserve the target definition. This
    # tests a physical invariance on one fixed training row, without refitting.
    doubled = counts[:1] * 2.0
    scaled = namespace["logcpm"](doubled) - np.column_stack([np.ones(1), namespace["conc_block"](doubled)]) @ coef
    scale_difference = float(np.max(np.abs(scaled[0] - reconstructed[0])))
    require(scale_difference <= ABS_TOL, "Target transform changes under uniform library scaling")
    normal_equations = design.T @ reconstructed
    relative_orthogonality = float(np.linalg.norm(normal_equations) / max(np.linalg.norm(design) * np.linalg.norm(logcpm), 1.0))
    require(relative_orthogonality <= 1e-12, "OLS residual orthogonality failed")
    args.output.mkdir(parents=True, exist_ok=False)
    np.save(args.output / "B99_concentration_coefficients.npy", coef, allow_pickle=False)
    np.save(args.output / "training_residual_mean.npy", mean_release, allow_pickle=False)
    for name, rows in (("participant_axis.tsv", participant_rows), ("h3k27ac_feature_axis.tsv", region_rows)):
        with (args.output / name).open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(rows)
    protocol = dict(schema_version="training-h3-target-transform-v1", descriptor_order=DESCRIPTORS,
                    target_formula="log2(1+1e6*C/C.sum(axis=1)) - column_stack([ones,conc_block(C)]) @ B99",
                    reference_n=99, n_regions=96460, coefficient_shape=[5, 96460], coefficient_dtype="float64",
                    entropy_base="natural", iqr_definition="numpy percentile75-minus25 of log2 positive counts",
                    top5_region_count=int(round(.05*96460)), region_semantics="one-based inclusive source keys; counted intervals [start-1,end)",
                    absolute_reproduction_tolerance_log2cpm=ABS_TOL,
                    tolerance_reason="float64 source/coefficients; permits small platform OLS roundoff, not float32 target approximation",
                    prediction_model_retrained=False, receiving_values_read=False, stochastic_operations=False)
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    for name, source in sources.items():
        (args.output / (name + "_source.txt")).write_text(source + "\n")
    report = dict(protocol=protocol, source_and_input_sha256=hashes, script_sha256=sha(Path(__file__)),
                  launcher_sha256=sha(Path(__file__).with_name("run_export_training_h3_target_transform.sbatch")),
                  checks=dict(rank=int(rank), design_condition_number=float(singular[0]/singular[-1]),
                              source_reproduction_max_absolute_difference=source_difference,
                              released_target_mean_max_absolute_difference=mean_difference,
                              uniform_scale_max_absolute_difference=scale_difference,
                              relative_residual_orthogonality=relative_orthogonality),
                  environment=dict(python=sys.version, numpy=np.__version__, platform=platform.platform(),
                                   slurm_job_id=os.environ["SLURM_JOB_ID"], threads={k:os.environ.get(k) for k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS")}),
                  output_sha256={p.name:sha(p) for p in sorted(args.output.iterdir()) if p.is_file()})
    (args.output / "receipt.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(dict(reference_n=99, n_regions=96460, coef_shape=list(coef.shape), **report["checks"]), sort_keys=True))


if __name__ == "__main__":
    main()
