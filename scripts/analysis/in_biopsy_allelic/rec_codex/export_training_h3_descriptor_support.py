"""Export training-only geometry for interpreting transported H3 residuals.

No prediction fit, receiving values, target recalibration or exclusion rule.
The saved geometry describes concentration descriptors used by the fixed B99
target definition. It is not a biological confidence region.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import platform
import sys

import numpy as np

from export_training_h3_target_transform import FIX, DESCRIPTORS, require, sha

REC = Path(__file__).resolve().parents[4] / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
TARGET = REC / "training_h3_target_transform_21997855"
INPUTS = {
    FIX / "h3k27ac_counts.npy": "9a4fa736a6dd5de0cd44e51428f2bf19e8dae9adb7a279e5b6c8f1c14103cd22",
    FIX / "participant_axis.tsv": "7e5be1a16df505e5daad451398e8b1a4b3cbbc48a32cd1381113e1ae4b7347bc",
    FIX / "h3k27ac_feature_axis.tsv": "cbe35aeb188e531283a1bcd636ba14de4cfbf215dfc544b99e3ebb97f9831986",
    TARGET / "conc_block_source.txt": "07740f7f5fa6660358926a64ca856a51ddc8a9aa247895f8e38d11ed81b446ed",
    TARGET / "B99_concentration_coefficients.npy": "8d34dafc551fb8b395c0d8b954707abfcc24097eb2aeb55479d99cf1950720a7",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Development molecular reads require compute job")
    require(not args.output.exists(), "Refusing output overwrite")
    hashes = {str(p): sha(p) for p in INPUTS}
    require(all(hashes[str(p)] == expected for p, expected in INPUTS.items()), "Frozen input changed")
    with (FIX / "participant_axis.tsv").open(newline="") as handle:
        participants = list(csv.DictReader(handle, delimiter="\t"))
    ids = [r["participant_id"] for r in participants]
    require(len(ids) == len(set(ids)) == 99, "Training biological roster changed")
    native = np.load(FIX / "h3k27ac_counts.npy", mmap_mode="r", allow_pickle=False)
    require(native.shape == (99, 96460) and native.dtype == np.uint32, "Training count axis changed")
    counts = np.asarray(native, dtype=np.float64)
    require(np.isfinite(counts).all() and (counts >= 0).all() and (counts.sum(1) > 0).all(), "Invalid counts")
    namespace = {"np": np}
    exec(compile((TARGET / "conc_block_source.txt").read_text(), "frozen-conc-block", "exec"), namespace)
    descriptors = namespace["conc_block"](counts)
    require(descriptors.shape == (99, 4) and np.isfinite(descriptors).all(), "Invalid descriptors")
    mean = descriptors.mean(0)
    scale = descriptors.std(0, ddof=0)
    require(np.isfinite(scale).all() and (scale > 0).all(), "Constant descriptor")
    standardized = (descriptors - mean) / scale
    design = np.column_stack([np.ones(99), standardized])
    u, singular, vt = np.linalg.svd(design, full_matrices=False)
    require(np.linalg.matrix_rank(design) == 5, "Descriptor geometry is rank deficient")
    # New-person leverage is sum(((z_new @ V) / singular_values)**2).
    # This stable factorization does not alter the raw-basis B99 coefficients.
    inverse_factor = vt.T / singular
    leverage = np.sum((design @ inverse_factor) ** 2, axis=1)
    qr, _ = np.linalg.qr(design, mode="reduced")
    qr_error = float(np.max(np.abs(leverage - np.sum(qr * qr, axis=1))))
    raw_design = np.column_stack([np.ones(99), descriptors])
    raw_qr, _ = np.linalg.qr(raw_design, mode="reduced")
    raw_error = float(np.max(np.abs(leverage - np.sum(raw_qr * raw_qr, axis=1))))
    require(qr_error <= 1e-12 and raw_error <= 1e-10, "Leverage depends on basis or differs from QR")
    require(abs(float(leverage.sum()) - 5.0) <= 1e-10, "Projection trace differs from rank")
    require((leverage >= 0).all() and (leverage <= 1 + 1e-12).all(), "Invalid training projection")
    scaled_descriptors = namespace["conc_block"](2.0 * counts[:1])
    scaled_design = np.r_[1.0, ((scaled_descriptors[0] - mean) / scale)]
    scale_leverage_error = abs(float(np.sum((scaled_design @ inverse_factor) ** 2)) - float(leverage[0]))
    require(scale_leverage_error <= 1e-10, "Leverage changes with uniform library scaling")
    args.output.mkdir(parents=True, exist_ok=False)
    np.savez(args.output / "training_descriptor_geometry.npz", descriptor_mean=mean,
             descriptor_scale=scale, descriptor_min=descriptors.min(0),
             descriptor_max=descriptors.max(0), inverse_design_factor=inverse_factor,
             training_leverage=leverage, standardized_design_singular_values=singular)
    with (args.output / "training_descriptors.tsv").open("x", newline="") as handle:
        fields = ["participant_id", *DESCRIPTORS[1:], "concentration_design_leverage"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for person, values, h in zip(ids, descriptors, leverage):
            writer.writerow(dict(participant_id=person, **dict(zip(DESCRIPTORS[1:], values)),
                                 concentration_design_leverage=float(h)))
    protocol = dict(schema_version="training-h3-descriptor-support-v1", biological_n=99,
                    n_regions=96460, descriptor_order=list(DESCRIPTORS[1:]),
                    descriptor_scaling="training mean and population standard deviation (ddof=0)",
                    new_person_leverage="sum((column_stack([1,(descriptors-mean)/scale]) @ inverse_design_factor)**2)",
                    interpretation="marginal ranges and leverage describe training descriptor geometry, not biological support probabilities",
                    prospective_receiving_report="all descriptor values, outside-training-range flags, leverage and ratio to maximum training leverage; retain every prespecified donor in primary error contrast",
                    exclusions_selected=False, receiving_coefficients_fitted=False,
                    prediction_model_retrained=False, receiving_values_read=False,
                    target_coefficients_changed=False, bootstrap_conditions_on_fixed_target_coefficients=True,
                    stochastic_operations=False)
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    summary = dict(protocol=protocol, source_and_input_sha256=hashes,
                   script_sha256=sha(Path(__file__)),
                   imported_helper_sha256=sha(Path(__file__).with_name("export_training_h3_target_transform.py")),
                   launcher_sha256=sha(Path(__file__).with_name("run_export_training_h3_descriptor_support.sbatch")),
                   checks=dict(qr_leverage_max_difference=qr_error, raw_basis_leverage_max_difference=raw_error,
                               uniform_scale_leverage_difference=scale_leverage_error,
                               training_leverage_sum=float(leverage.sum()), rank=5),
                   descriptor_min=dict(zip(DESCRIPTORS[1:], descriptors.min(0))),
                   descriptor_max=dict(zip(DESCRIPTORS[1:], descriptors.max(0))),
                   max_training_leverage=float(leverage.max()),
                   standardized_design_condition=float(singular[0] / singular[-1]),
                   environment=dict(python=sys.version, numpy=np.__version__, platform=platform.platform(),
                                    slurm_job_id=os.environ["SLURM_JOB_ID"]),
                   output_sha256={p.name: sha(p) for p in sorted(args.output.iterdir()) if p.is_file()})
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["checks"], sort_keys=True))


if __name__ == "__main__":
    main()
