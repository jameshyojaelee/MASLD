"""Audit fixed short-reference features in permitted development RNA/model.

No receiving data, prediction, refit, gene exclusion, or zero filling. Direct
log-CPM slopes do not remove any feature's role in the common denominator.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np

from export_training_h3_target_transform import FIX, RELEASE, require, sha, read_axis

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
SHORT = REC / "ensembl98_salmon_metadata_21998488/modeled_genes_no_unclipped_k31_transcript.tsv"
INPUTS = {
    FIX / "rna_values.npy": "812d5fce6e7cc88f6939ae493db7b7e7470057c21cf3b219abc39f587ff64ddc",
    FIX / "rna_feature_axis.tsv": "baf983ad18a893b7e5db04364818073facc68b47ae25f63070373481e3b217dd",
    FIX / "participant_axis.tsv": "7e5be1a16df505e5daad451398e8b1a4b3cbbc48a32cd1381113e1ae4b7347bc",
    RELEASE / "weights/chromatin_state_v1_2.npz": "72f9eaba359f7c4c61cac2e00da2d5539b5d3ec868dc33b1bacc2ef0105694ec",
    SHORT: "4f348047d98367f9e041716dbe9a7a8952adcf54bba78c0dc97f55477151f7bf",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Development RNA/model audit requires compute")
    require(not args.output.exists(), "Refusing output overwrite")
    for path, expected in INPUTS.items():
        require(sha(path) == expected, "Frozen input changed: " + str(path))
    axis = read_axis(FIX / "rna_feature_axis.tsv")
    people = read_axis(FIX / "participant_axis.tsv")
    short = read_axis(SHORT)
    genes = np.array([r["stable_gene_id"] for r in axis])
    require(len(genes) == len(set(genes)) == 42163, "Modeled gene axis changed")
    require(len(people) == len({r["participant_id"] for r in people}) == 99,
            "Development biological roster changed")
    require(len(short) == 40, "Short-gene population changed")
    indices = np.array([int(r["rna_feature_index"]) for r in short])
    require(np.array_equal(genes[indices], [r["stable_gene_id"] for r in short]),
            "Short-gene axis join differs")
    counts = np.load(FIX / "rna_values.npy", mmap_mode="r", allow_pickle=False)
    require(counts.shape == (99, 42163) and counts.dtype == np.float64, "RNA shape/dtype changed")
    require(np.isfinite(counts).all() and (counts >= 0).all(), "Invalid RNA")
    library = counts.sum(1)
    require((library > 0).all(), "Empty development RNA library")
    values = np.asarray(counts[:, indices])
    mass = values / library[:, None]
    with np.load(RELEASE / "weights/chromatin_state_v1_2.npz", allow_pickle=True) as weights:
        require(np.array_equal(weights["gene_axis"].astype(str), genes), "Packaged RNA axis changed")
        require(weights["prof_region_key"].shape == (96460,), "Packaged H3 axis changed")
        xmean = np.asarray(weights["prof_x_mean"], float)
        xs = np.asarray(weights["prof_rrr_xs"], float)
        xsd = np.asarray(weights["prof_x_sd"], float)
        require(xs.shape == xsd.shape == (42163,) and (xs > 0).all() and (xsd > 0).all(),
                "Invalid fixed RNA scales")
        source_mean = np.log2(1 + 1e6 * counts / library[:, None]).mean(0)
        mean_error = float(np.max(np.abs(source_mean - xmean)))
        require(mean_error <= 1e-10, "Development RNA differs from released training mean")
        # At fixed library denominator, these are direct slopes per log2-CPM.
        # Construct only 40 x 96,460 rather than the whole gene x region matrix.
        global_slopes = ((np.asarray(weights["prof_rrr_VtDv"], float)[indices]
                          / xs[indices, None]) @ np.asarray(weights["prof_rrr_WVr"], float))
        global_slopes = global_slopes @ np.asarray(weights["prof_rrr_Vr"], float).T
        require(global_slopes.shape == (40, 96460) and np.isfinite(global_slopes).all(),
                "Invalid effective global slopes")
        cis_index = np.asarray(weights["prof_Ocis_cis_idx"])
        cis_coef = np.asarray(weights["prof_Ocis_coef_cis"], float)
        global_head = np.asarray(weights["prof_Ocis_coef_glob"], float)
        require(not np.any(global_head), "Ocis unexpectedly depends on state intermediates")
        require(cis_index.shape == cis_coef.shape and cis_index.shape[0] == 96460,
                "Invalid local coefficient geometry")
        require(cis_index.dtype.kind in "iu" and ((cis_index >= -1) & (cis_index < 42163)).all(),
                "Invalid local feature indices")
        require(np.isfinite(cis_coef).all(), "Nonfinite local coefficients")
        primary_slopes = global_slopes.copy()
        records = []
        for j, (row, index) in enumerate(zip(short, indices)):
            regions, slots = np.nonzero(cis_index == index)
            slopes = cis_coef[regions, slots] / xsd[index]
            np.add.at(primary_slopes[j], regions, slopes)
            records.append(dict(stable_gene_id=row["stable_gene_id"], rna_feature_index=int(index),
                                development_donors_positive=int(np.count_nonzero(values[:, j])),
                                development_count_sum=float(values[:, j].sum()),
                                development_count_max=float(values[:, j].max()),
                                development_mean_modeled_rna_mass=float(mass[:, j].mean()),
                                development_max_modeled_rna_mass=float(mass[:, j].max()),
                                fixed_denominator_member=True,
                                global_direct_slope_nonzero=bool(np.any(global_slopes[j])),
                                global_direct_slope_max_abs=float(np.abs(global_slopes[j]).max()),
                                local_nonzero_region_count=int(len(set(regions[slopes != 0]))),
                                local_direct_slope_max_abs=float(np.abs(slopes).max()) if len(slopes) else 0.0,
                                primary_direct_slope_nonzero=bool(np.any(primary_slopes[j])),
                                primary_direct_slope_max_abs=float(np.abs(primary_slopes[j]).max())))
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "short_gene_development_activity.tsv").open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    summary = dict(biological_n=99, modeled_genes=42163, short_modeled_genes=40,
                   short_genes_positive_in_development=sum(r["development_donors_positive"] > 0 for r in records),
                   short_genes_nonzero_global_direct_slopes=sum(r["global_direct_slope_nonzero"] for r in records),
                   short_genes_nonzero_local_direct_slopes=sum(r["local_nonzero_region_count"] > 0 for r in records),
                   short_genes_nonzero_primary_direct_slopes=sum(r["primary_direct_slope_nonzero"] for r in records),
                   short_set_rna_mass_min=float(mass.sum(1).min()),
                   short_set_rna_mass_mean=float(mass.sum(1).mean()),
                   short_set_rna_mass_max=float(mass.sum(1).max()),
                   training_mean_max_difference=mean_error,
                   direct_slope_units="prediction residual log2-CPM per input log2-CPM, holding modeled-gene library denominator fixed",
                   denominator_coupling="every modeled gene contributes to the sum before log-CPM; direct zero slopes do not establish safe omission",
                   zero_interpretation="development zeros do not establish biological absence, receiving zeros, or complete RNA observability",
                   receiving_data_read=False, model_fitted=False, predictions_computed=False,
                   exclusions_or_zero_fill=False, contract_relaxed=False, inference_or_multiple_tests=False,
                   input_sha256={str(p): h for p, h in INPUTS.items()},
                   script_sha256=sha(Path(__file__)),
                   helper_sha256=sha(Path(__file__).with_name("export_training_h3_target_transform.py")),
                   launcher_sha256=sha(Path(__file__).with_name("run_audit_short_rna_model_features.sbatch")),
                   environment=dict(python=sys.version, numpy=np.__version__, slurm_job_id=os.environ["SLURM_JOB_ID"]),
                   output_sha256={p.name: sha(p) for p in args.output.iterdir() if p.is_file()})
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k.startswith("short_")}, sort_keys=True))


if __name__ == "__main__":
    main()
