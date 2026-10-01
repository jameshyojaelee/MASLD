"""Apply the fixed B99 H3 residual definition without fitting receiving data.

The command performs synthetic checks only. The callable requires already
admitted integer single-read counts on all seventeen donors and 96,460 regions.
It cannot establish donor identity, count provenance or assay comparability.
"""
import ast
import json
import os
from pathlib import Path
import sys

import numpy as np

from codex_paired_profile_evaluation import DONORS, REGIONS
from export_training_h3_target_transform import DESCRIPTORS, FIX, read_axis, require, sha

REC = Path(__file__).resolve().parents[4] / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
TARGET = REC / "training_h3_target_transform_21997855"
SUPPORT = REC / "training_h3_descriptor_support_21998378"
GUARDS = {
    FIX / "h3k27ac_feature_axis.tsv": "cbe35aeb188e531283a1bcd636ba14de4cfbf215dfc544b99e3ebb97f9831986",
    TARGET / "B99_concentration_coefficients.npy": "8d34dafc551fb8b395c0d8b954707abfcc24097eb2aeb55479d99cf1950720a7",
    TARGET / "training_residual_mean.npy": "a16179d29b2a45ccb374987796ab75393528c54fc38f255b30faa6af255b7fa3",
    TARGET / "conc_block_source.txt": "07740f7f5fa6660358926a64ca856a51ddc8a9aa247895f8e38d11ed81b446ed",
    TARGET / "protocol.json": "3f5d236523dd5155e079f5e002b3eab13c197c81e94500534711c6824a42979a",
    SUPPORT / "training_descriptor_geometry.npz": "43a06ca9d1591961664fb78c477560a3cfd91c0c65bf08a7fab6a7f4b17da5e9",
    SUPPORT / "protocol.json": "f0d58e8d1567aa806886b78e8cfcab3040a0c5a501a80658ae176b05d5d81d9e",
}


def load_fixed_definition():
    require(os.environ.get("SLURM_JOB_ID"), "Load fixed numerical artifacts on compute")
    for path, expected in GUARDS.items():
        require(sha(path) == expected, "Frozen target/support artifact changed: " + str(path))
    protocol = json.loads((TARGET / "protocol.json").read_text())
    require(tuple(protocol["descriptor_order"]) == DESCRIPTORS
            and protocol["reference_n"] == 99 and protocol["n_regions"] == REGIONS,
            "Target descriptor convention changed")
    keys = tuple(row["opaque_source_feature_key"] for row in read_axis(FIX / "h3k27ac_feature_axis.tsv"))
    require(len(keys) == len(set(keys)) == REGIONS, "Frozen region identity differs")
    coef = np.load(TARGET / "B99_concentration_coefficients.npy", allow_pickle=False)
    mean = np.load(TARGET / "training_residual_mean.npy", allow_pickle=False)
    require(coef.shape == (5, REGIONS) and mean.shape == (REGIONS,)
            and coef.dtype == mean.dtype == np.float64
            and np.isfinite(coef).all() and np.isfinite(mean).all(), "Invalid fixed target coefficients")
    with np.load(SUPPORT / "training_descriptor_geometry.npz", allow_pickle=False) as archive:
        geometry = {name: archive[name].copy() for name in archive.files}
    require(all(np.isfinite(value).all() for value in geometry.values()), "Nonfinite training geometry")
    require(all(geometry[name].shape == (4,) for name in
                ("descriptor_mean", "descriptor_scale", "descriptor_min", "descriptor_max"))
            and geometry["inverse_design_factor"].shape == (5, 5)
            and geometry["training_leverage"].shape == (99,)
            and (geometry["descriptor_scale"] > 0).all(), "Invalid training geometry axes")
    source = (TARGET / "conc_block_source.txt").read_text()
    tree = ast.parse(source)
    require(len(tree.body) == 1 and isinstance(tree.body[0], ast.FunctionDef)
            and tree.body[0].name == "conc_block", "Frozen descriptor source is not one function")
    namespace = {"np": np}
    exec(compile(tree, str(TARGET / "conc_block_source.txt"), "exec"), namespace)
    return dict(region_keys=keys, coefficients=coef, training_mean=mean,
                geometry=geometry, concentration=namespace["conc_block"])


def apply_fixed_target(counts, donor_ids, region_keys):
    """Return residual labels and descriptor diagnostics; retain every donor.

    Identity here means supplied axis equality. Actual sample/count provenance
    must be established by the acquisition and counting procedure beforehand.
    Diagnostic range/leverage flags never select or omit receiving donors.
    """
    require(os.environ.get("SLURM_JOB_ID"), "Numerical target transformation requires compute")
    require(tuple(donor_ids) == DONORS, "Require all seventeen fixed donors in order")
    definition = load_fixed_definition()
    require(tuple(region_keys) == definition["region_keys"], "Require the complete fixed region identity/order")
    counts = np.asarray(counts)
    require(counts.shape == (len(DONORS), REGIONS)
            and np.issubdtype(counts.dtype, np.integer), "Require integer single-read counts on the entire axis")
    require((counts >= 0).all() and int(counts.max()) <= 2**53, "Counts outside exact nonnegative float64 integer range")
    # Python integers avoid uint64 overflow and floating-point rounding at
    # the boundary. All partial nonnegative sums then remain exactly representable.
    exact_libraries = [sum(map(int, row)) for row in counts]
    require(all(0 < value <= 2**53 for value in exact_libraries),
            "Every donor requires a positive region-count sum within the exact float64 integer range")
    raw = counts.astype(np.float64)
    libraries = np.asarray(exact_libraries, dtype=np.float64)
    descriptors = definition["concentration"](raw)
    require(descriptors.shape == (len(DONORS), 4) and np.isfinite(descriptors).all(), "Invalid receiving concentration descriptors")
    design = np.column_stack([np.ones(len(DONORS)), descriptors])
    residual = np.log2(1.0 + 1e6 * raw / libraries[:, None]) - design @ definition["coefficients"]
    require(residual.dtype == np.float64 and np.isfinite(residual).all(), "Nonfinite fixed residual labels")
    geometry = definition["geometry"]
    standardized = np.column_stack([np.ones(len(DONORS)),
                                   (descriptors - geometry["descriptor_mean"]) / geometry["descriptor_scale"]])
    leverage = np.sum((standardized @ geometry["inverse_design_factor"])**2, axis=1)
    maximum_training_leverage = float(geometry["training_leverage"].max())
    require(maximum_training_leverage > 0 and np.isfinite(leverage).all(), "Invalid descriptor leverage")
    outside = (descriptors < geometry["descriptor_min"]) | (descriptors > geometry["descriptor_max"])
    diagnostics = []
    for i, donor in enumerate(DONORS):
        diagnostics.append(dict(donor_id=donor, region_count_sum=float(libraries[i]),
                                descriptors=dict(zip(DESCRIPTORS[1:], descriptors[i].tolist())),
                                outside_training_marginal_range=dict(zip(DESCRIPTORS[1:], outside[i].tolist())),
                                any_descriptor_outside_training_range=bool(outside[i].any()),
                                concentration_design_leverage=float(leverage[i]),
                                leverage_over_maximum_training=float(leverage[i] / maximum_training_leverage)))
    return residual, diagnostics


def check_invariants():
    require(os.environ.get("SLURM_JOB_ID"), "Synthetic numerical checks require compute")
    definition = load_fixed_definition()
    columns = np.arange(REGIONS, dtype=np.uint32)
    counts = np.stack([1 + (columns % 97) * (i % 5 + 1) for i in range(len(DONORS))])
    counts[:, ::11] = 0
    counts[0] = 1
    counts[1] = 1
    counts[1, 0] = 1_000_000
    target, diagnostics = apply_fixed_target(counts, DONORS, definition["region_keys"])
    scaled, scaled_diagnostics = apply_fixed_target(2 * counts, DONORS, definition["region_keys"])
    scale_error = float(np.max(np.abs(target - scaled)))
    leverage_error = max(abs(a["concentration_design_leverage"] - b["concentration_design_leverage"])
                         for a, b in zip(diagnostics, scaled_diagnostics))
    require(scale_error <= 1e-9 and leverage_error <= 1e-10, "Uniform library scaling changes the fixed target or descriptor geometry")
    # Independently derived uniform-count row: Gini=0, entropy=ln(R),
    # positive-log-count IQR=0, and top5 fraction=round(.05R)/R.
    uniform = np.array([1., 0., np.log(REGIONS), 0., round(.05 * REGIONS) / REGIONS])
    expected = np.log2(1. + 1e6 / REGIONS) - uniform @ definition["coefficients"]
    uniform_error = float(np.max(np.abs(target[0] - expected)))
    require(uniform_error <= 1e-9, "Hand-derived uniform-count residual differs")
    require(len(diagnostics) == len(DONORS)
            and diagnostics[0]["any_descriptor_outside_training_range"]
            and diagnostics[1]["any_descriptor_outside_training_range"], "Range diagnostics dropped extreme synthetic donors")
    rounded_overflow = counts.astype(np.uint64)
    rounded_overflow[0] = 0
    rounded_overflow[0, :2] = [2**53, 1]
    invalid = [(counts, DONORS[::-1], definition["region_keys"]),
               (counts, DONORS, definition["region_keys"][::-1]),
               (counts, DONORS, definition["region_keys"][:-1]),
               (counts[:-1], DONORS, definition["region_keys"]),
               (np.zeros_like(counts), DONORS, definition["region_keys"]),
               (counts.astype(float) + .5, DONORS, definition["region_keys"]),
               (-counts.astype(np.int64), DONORS, definition["region_keys"]),
               (rounded_overflow, DONORS, definition["region_keys"])]
    for values, donors, regions in invalid:
        try:
            apply_fixed_target(values, donors, regions)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid count/axis input was accepted")
    for path, expected_sha in GUARDS.items():
        require(sha(path) == expected_sha, "Frozen artifact changed during synthetic checks")
    output = REC / ("fixed_h3_target_checks_" + os.environ["SLURM_JOB_ID"])
    output.mkdir(exist_ok=False)
    summary = dict(uniform_scale_target_max_difference=scale_error,
                   uniform_scale_leverage_max_difference=leverage_error,
                   uniform_count_hand_derivation_max_difference=uniform_error,
                   invalid_input_refusals=len(invalid), retained_synthetic_donors=len(DONORS),
                   synthetic_descriptor_diagnostics=diagnostics,
                   input_sha256={str(p): v for p, v in GUARDS.items()},
                   script_sha256=sha(Path(__file__)),
                   launcher_sha256=sha(Path(__file__).with_name("run_fixed_h3_target_checks.sbatch")),
                   protocol=dict(target="fixed B99 concentration-residualized log2 modeled-region CPM",
                                 full_regions=REGIONS, receiving_donors=list(DONORS),
                                 count_unit="admitted primary nonsupplementary single reads, MAPQ>=10, unstranded, featureCounts without -O/-M; no unique-mapping guarantee",
                                 fitted_coefficients_changed=False, donors_selected_by_diagnostics=False,
                                 interpretation="descriptor geometry is not a biological confidence interval"),
                   receiving_values_read=False, prediction_accuracy_evaluated=False,
                   biological_validation=False, stochastic_operations=False,
                   environment=dict(python=sys.version, numpy=np.__version__, slurm_job_id=os.environ["SLURM_JOB_ID"]))
    (output / "executed_source.py").write_bytes(Path(__file__).read_bytes())
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("uniform_scale_target_max_difference", "uniform_scale_leverage_max_difference", "uniform_count_hand_derivation_max_difference", "invalid_input_refusals", "retained_synthetic_donors")}))


if __name__ == "__main__":
    check_invariants()
