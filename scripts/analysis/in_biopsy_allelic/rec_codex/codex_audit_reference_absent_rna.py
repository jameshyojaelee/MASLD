"""Development RNA relevance of genes absent from official Ensembl98 FASTAs.

No receiving data, predictor fitting, sequence filtering, or zero filling.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np

from export_training_h3_target_transform import FIX, require, sha, read_axis

ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / "Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z"
MISSING = REC / "ensembl98_reference_census_21998377/missing_modeled_union_genes.tsv"
INPUTS = {
    FIX / "rna_values.npy": "812d5fce6e7cc88f6939ae493db7b7e7470057c21cf3b219abc39f587ff64ddc",
    FIX / "rna_feature_axis.tsv": "baf983ad18a893b7e5db04364818073facc68b47ae25f63070373481e3b217dd",
    FIX / "participant_axis.tsv": "7e5be1a16df505e5daad451398e8b1a4b3cbbc48a32cd1381113e1ae4b7347bc",
    MISSING: "a922ce25bd465d906f27730b02a4a475d3ee56730c8a2c96fb64f72eef52cadd",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get("SLURM_JOB_ID"), "Development molecular audit requires compute")
    require(not args.output.exists(), "Refusing output overwrite")
    for p, expected in INPUTS.items():
        require(sha(p) == expected, "Frozen input changed: " + str(p))
    axis = read_axis(FIX / "rna_feature_axis.tsv")
    people = read_axis(FIX / "participant_axis.tsv")
    missing = read_axis(MISSING)
    require(len(people) == len({r["participant_id"] for r in people}) == 99, "Donor roster changed")
    require(len(axis) == len({r["stable_gene_id"] for r in axis}) == 42163, "RNA axis changed")
    require(len(missing) == 1010, "Reference absence population changed")
    indices = np.array([int(r["rna_feature_index"]) for r in missing])
    require(len(set(indices)) == 1010 and (indices >= 0).all() and (indices < 42163).all(), "Invalid absent-gene indices")
    for r, i in zip(missing, indices):
        require(axis[i]["stable_gene_id"] == r["stable_gene_id"] and
                r["modeled_in_union"] == "False" and r["modeled_in_gtf"] == "True",
                "Reference absence/gene identity differs")
    counts = np.load(FIX / "rna_values.npy", mmap_mode="r", allow_pickle=False)
    require(counts.shape == (99, 42163) and counts.dtype == np.float64, "RNA shape/dtype differs")
    require(np.isfinite(counts).all() and (counts >= 0).all(), "Invalid RNA estimates")
    library = counts.sum(1)
    require(np.isfinite(library).all() and (library > 0).all(), "Invalid modeled-gene library total")
    values = np.asarray(counts[:, indices])
    mass = values.sum(1) / library
    rows = [dict(stable_gene_id=r["stable_gene_id"], rna_feature_index=int(i),
                 donors_with_positive_training_estimate=int(np.count_nonzero(values[:, j])),
                 training_estimate_sum=float(values[:, j].sum()),
                 training_estimate_max=float(values[:, j].max()))
            for j, (r, i) in enumerate(zip(missing, indices))]
    args.output.mkdir(parents=True, exist_ok=False)
    path = args.output / "reference_absent_gene_training_estimates.tsv"
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary = dict(biological_n=99, reference_absent_modeled_genes=1010,
                   genes_positive_in_training=sum(r["donors_with_positive_training_estimate"] > 0 for r in rows),
                   genes_zero_in_all_training_donors=sum(r["donors_with_positive_training_estimate"] == 0 for r in rows),
                   modeled_library_fraction_of_absent_set=dict(min=float(mass.min()), median=float(np.median(mass)), max=float(mass.max())),
                   estimate_semantics="native fractional source gene estimates; zero does not establish biological absence or receiving zero",
                   reference_semantics="absent from official cDNA/ncRNA FASTAs but present in matching GTF; reconstructed reference supplies native annotated transcripts",
                   receiving_data_read=False, model_fitted=False, predictions_computed=False,
                   genes_removed_or_filled=False, private_pipeline_equivalence_proved=False,
                   input_sha256={str(p): v for p, v in INPUTS.items()},
                   script_sha256=sha(Path(__file__)),
                   helper_sha256=sha(Path(__file__).with_name("export_training_h3_target_transform.py")),
                   launcher_sha256=sha(Path(__file__).with_name("run_codex_audit_reference_absent_rna.sbatch")),
                   output_sha256={path.name: sha(path)},
                   environment=dict(python=sys.version, numpy=np.__version__, slurm_job_id=os.environ["SLURM_JOB_ID"]))
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({k: summary[k] for k in ("genes_positive_in_training", "genes_zero_in_all_training_donors", "modeled_library_fraction_of_absent_set")}, sort_keys=True))


if __name__ == "__main__":
    main()
