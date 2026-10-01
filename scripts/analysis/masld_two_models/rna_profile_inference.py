#!/usr/bin/env python3
"""Candidate RNA-to-profile interface with an optional fixed input-abundance reference.

Inputs are a genes-by-samples TSV of counts or the continuous deposited RNA
estimates. TPM, log-expression and normalized expression units are unsupported.
The reference-mass transformation was examined only with the bounded rank-12,
1,000-region development recipe. Its use with the full released weights here
is an inference extension without held-participant accuracy evaluation.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import platform
import sys
import time
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
RELEASE = ROOT / "Analysis/MASLD_Model_Benchmark/release/masld-liver-chromatin-state-v1.2"
WEIGHTS = RELEASE / "weights/chromatin_state_v1_2.npz"
FIX = ROOT / "Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture"
CANDIDATE = ROOT / "GWAS/finemapping/results/alphagenome_campaign/two-models-20260922T203900EDT"
DEFAULT_REFERENCE = CANDIDATE / "rna_profile_reference_candidate"
FORM = "rrr_offset_cis"
REFERENCE_MASS_DEFINITION = (
    "Share of summed continuous input abundance among the 42,163 modeled genes, "
    "averaged with equal weight across training participants. This is not physical "
    "RNA mass, a molecule fraction, or a whole-transcriptome fraction."
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@lru_cache(maxsize=1)
def released_scorer():
    spec = importlib.util.spec_from_file_location("candidate_released_profile_scorer", RELEASE / "score.py")
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    return scorer


def source_column(path, column):
    with Path(path).open(newline="") as stream:
        return np.array([row[column] for row in csv.DictReader(stream, delimiter="\t")], dtype=str)


def export_reference(destination):
    """Freeze the release-training reference without reading any chromatin label."""
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Export the training reference within a compute allocation")
    if destination.exists():
        raise FileExistsError(destination)
    scorer = released_scorer()
    participant_file = FIX / "molecular/participant_axis.tsv"
    gene_file = FIX / "molecular/rna_feature_axis.tsv"
    expression_file = FIX / "molecular/rna_values.npy"
    participants = source_column(participant_file, "participant_id")
    genes = scorer.normalise_ids(source_column(gene_file, "stable_gene_id"))
    values = np.asarray(np.load(expression_file, mmap_mode="r"), float)
    with np.load(WEIGHTS, allow_pickle=True) as weights:
        axis = weights["gene_axis"].astype(str)
        training_mean = np.asarray(weights["prof_x_mean"], float)
        training_n = int(weights["n_train"])
    if training_n != 99 or len(participants) != 99 or len(set(participants)) != 99:
        raise ValueError("Release training participant census differs")
    if values.shape != (len(participants), len(genes)) or len(set(genes)) != len(genes):
        raise ValueError("Training expression gene axis is ambiguous")
    if set(genes) != set(axis):
        raise ValueError("Training and released gene universes differ")
    positions = {gene: i for i, gene in enumerate(genes)}
    values = values[:, [positions[gene] for gene in axis]]
    if not np.isfinite(values).all() or np.any(values < 0) or np.any(values.sum(1) <= 0):
        raise ValueError("Invalid continuous training expression")
    proportions = values / values.sum(1)[:, None]
    expected_mass = proportions.mean(0)
    mean_difference = float(np.max(np.abs(np.log2(proportions * 1e6 + 1).mean(0) - training_mean)))
    if mean_difference > 1e-10:
        raise ValueError("Training RNA does not reconstruct the released expression means")
    destination.mkdir(parents=True)
    artifact = destination / "reference.npz"
    np.savez_compressed(artifact, gene_axis=axis, expected_gene_RNA_mass=expected_mass,
                        participant_id=participants)
    report = {
        "status": "candidate fixed RNA reference; no model refit or named-release change",
        "n_training_participants": len(participants), "n_genes": len(axis),
        "reference": "equal-participant mean of modeled-gene input-abundance proportions",
        "reference_mass_definition": REFERENCE_MASS_DEFINITION,
        "input_units": "continuous deposited estimates, never rounded",
        "training_log2cpm_mean_max_difference_from_release": mean_difference,
        "source_sha256": {str(path): sha256(path) for path in
                          (WEIGHTS, RELEASE / "score.py", expression_file, gene_file, participant_file)},
        "reference_npz_sha256": sha256(artifact),
        "scope": "Modeled-gene input-abundance proportions from the 99 release-training participants. Not a new evaluation source. Full-weight reference-mass prediction has no held-participant accuracy evaluation.",
    }
    (destination / "reference.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def load_reference(directory):
    artifact = directory / "reference.npz"
    report = json.loads((directory / "reference.json").read_text())
    if sha256(artifact) != report["reference_npz_sha256"]:
        raise ValueError("Candidate RNA reference contents changed")
    for source in (WEIGHTS, RELEASE / "score.py"):
        if sha256(source) != report["source_sha256"][str(source)]:
            raise ValueError("Release weights or scoring code differ from the fixed RNA reference")
    with np.load(artifact, allow_pickle=False) as source:
        axis = source["gene_axis"].astype(str)
        mass = np.asarray(source["expected_gene_RNA_mass"], float)
    if mass.shape != axis.shape or not np.isfinite(mass).all() or np.any(mass < 0) or not np.isclose(mass.sum(), 1, atol=1e-12):
        raise ValueError("Invalid reference RNA proportions")
    return axis, mass, report


def prepare_input(counts, genes, axis, expected_mass):
    scorer = released_scorer()
    counts = np.asarray(counts, float)
    if counts.ndim != 2 or counts.shape[0] != len(genes) or counts.shape[1] == 0:
        raise ValueError("Expected genes-by-samples RNA matrix")
    if not np.isfinite(counts).all() or np.any(counts < 0):
        raise ValueError("RNA values must be finite and nonnegative")
    normalized = scorer.normalise_ids(genes)
    counts, normalized, duplicates = scorer.collapse_duplicates(counts, normalized)
    positions = {gene: i for i, gene in enumerate(normalized)}
    hit = np.array([i for i, gene in enumerate(axis) if gene in positions], dtype=int)
    if len(hit) / len(axis) < scorer.MIN_COVERAGE:
        raise ValueError(f"Gene coverage below the inherited release input floor {scorer.MIN_COVERAGE}; this is not an accuracy guarantee")
    matched = np.array([positions[axis[i]] for i in hit], dtype=int)
    observed = counts[matched]
    library = observed.sum(0)
    if np.any(library <= 0):
        raise ValueError("Every sample requires positive observed RNA totals")
    if float(counts.max()) < scorer.LOG_LIKE_MAX:
        raise ValueError("Input appears log transformed; use counts or continuous deposited estimates")
    # Enforce exact equality to the released raw route for the complete axis.
    covered_mass = 1.0 if len(hit) == len(axis) else float(expected_mass[hit].sum())
    if not 0 < covered_mass <= 1:
        raise ValueError("Observed genes have invalid expected reference RNA mass")
    coverage = {
        "n_axis_genes": len(axis), "n_matched_genes": len(hit),
        "gene_coverage_fraction": len(hit) / len(axis),
        "expected_observed_reference_RNA_mass_fraction": covered_mass,
        "reference_mass_definition": REFERENCE_MASS_DEFINITION,
        "n_duplicate_stable_id_rows_summed": int(duplicates),
        "n_unmatched_input_genes": len(normalized) - len(hit),
        "matched_input_total_per_sample": library.tolist(),
        "input_floor": "Inherited release minimum gene fraction; no new threshold selected and no accuracy certification",
    }
    return observed, hit, covered_mass, coverage


def reference_mass_input(observed, hit, covered_mass, training_mean):
    library = observed.sum(0)
    cpm = observed / library[None, :] * 1e6
    if covered_mass != 1:
        cpm *= covered_mass
    transformed = np.tile(training_mean, (observed.shape[1], 1))
    transformed[:, hit] = np.log2(cpm + 1).T
    return transformed


def apply_reference_mass_weights(transformed, weights):
    """Apply only the verified released RRR plus local-gene residual form."""
    region = weights["prof_region_key"].astype(str)
    global_coefficients = np.asarray(weights["prof_Ocis_coef_glob"], float)
    if global_coefficients.shape != (len(region), 0):
        raise ValueError("The local head contains global/state contributions; this candidate cannot silently omit them")
    xmean = np.asarray(weights["prof_x_mean"], float)
    xsd = np.asarray(weights["prof_x_sd"], float)
    rrr_sd = np.asarray(weights["prof_rrr_xs"], float)
    if np.any(xsd <= 0) or np.any(rrr_sd <= 0):
        raise ValueError("Invalid released RNA scales")
    rrr = ((((transformed - np.asarray(weights["prof_rrr_xm"], float)) / rrr_sd)
            @ np.asarray(weights["prof_rrr_VtDv"], float))
           @ np.asarray(weights["prof_rrr_WVr"], float))
    prediction = rrr @ np.asarray(weights["prof_rrr_Vr"], float).T + np.asarray(weights["prof_rrr_ym"], float)
    local = np.tile(np.asarray(weights["prof_Ocis_intercept"], float), (len(transformed), 1))
    standardized = (transformed - xmean) / xsd
    indices = np.asarray(weights["prof_Ocis_cis_idx"], int)
    coefficient = np.asarray(weights["prof_Ocis_coef_cis"], float)
    if indices.shape != coefficient.shape or indices.shape[0] != len(region) or np.any(indices >= len(xmean)) or np.any(indices < -1):
        raise ValueError("Local-head gene indices or coefficients differ")
    for k in range(indices.shape[1]):
        present = indices[:, k] >= 0
        local[:, present] += standardized[:, indices[present, k]] * coefficient[present, k][None, :]
    prediction += local
    if not np.isfinite(prediction).all():
        raise ValueError("Nonfinite predicted profile")
    return prediction, region


def score_profiles(counts, genes, units, transport="raw", reference=DEFAULT_REFERENCE):
    if units not in ("counts", "deposited_estimates"):
        raise ValueError("Declare counts or deposited_estimates; TPM and normalized/log units are unsupported")
    if transport not in ("raw", "marginal", "reference_mass"):
        raise ValueError("Unknown input transformation")
    start = time.perf_counter()
    axis, expected_mass, reference_report = load_reference(Path(reference))
    observed, hit, covered_mass, coverage = prepare_input(counts, genes, axis, expected_mass)
    with np.load(WEIGHTS, allow_pickle=True) as weights:
        if not np.array_equal(axis, weights["gene_axis"].astype(str)):
            raise ValueError("Candidate reference and released gene order differ")
        if transport == "reference_mass":
            transformed = reference_mass_input(observed, hit, covered_mass,
                                               np.asarray(weights["prof_x_mean"], float))
            prediction, region = apply_reference_mass_weights(transformed, weights)
        else:
            released = released_scorer().profile(counts, genes, weights_path=str(WEIGHTS),
                                                 form=FORM, transport=transport)
            prediction, region = released["profile"], released["region_key"]
    report = {
        "status": "candidate research inference", "form": FORM, "transport": transport,
        "input_units_declared": units, "input_measurements_rounded": False,
        "output_units": "predicted concentration-residual H3K27ac log2 CPM in GSE267145 training units",
        "profile_is_measured_chromatin": False, "n_samples": prediction.shape[0],
        "n_modeled_regions": prediction.shape[1], "coverage": coverage,
        "individual_region_accuracy_claim": False, "mechanism_claim": False,
        "individual_predictive_intervals": None,
        "reference_mass_scope": "Adaptive development evidence used a fixed rank-12/1,000-region recipe in 99 source participants. The full released rank-20/96,460-region inference extension implemented here has no held-participant accuracy evaluation.",
        "reference_mass_definition": REFERENCE_MASS_DEFINITION,
        "reference_mass_assumption": "The observed gene set's share of summed input abundance over the modeled genes is compatible with the fixed training reference; gene coverage and expected abundance coverage are not individual accuracy guarantees.",
        "reference_npz_sha256": reference_report["reference_npz_sha256"],
        "weights_sha256": reference_report["source_sha256"][str(WEIGHTS)],
        "elapsed_seconds_including_reference_validation": time.perf_counter() - start,
        "environment": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform(),
                        "slurm_job_id": os.environ.get("SLURM_JOB_ID")},
    }
    return {"profile": prediction, "region_key": region, "report": report}


def main(args):
    if args.export_reference is not None:
        if args.counts is not None or args.out is not None:
            raise ValueError("Export the fixed reference separately from scoring")
        report = export_reference(args.export_reference)
        print(json.dumps({"reference": str(args.export_reference), "n_genes": report["n_genes"]}))
        return
    if args.counts is None or args.units is None or args.out is None:
        raise ValueError("Scoring requires --counts, --units and --out")
    if args.out.exists():
        raise FileExistsError(args.out)
    scorer = released_scorer()
    counts, genes, samples = scorer.read_counts(str(args.counts))
    if len(set(samples)) != len(samples) or any(not name for name in samples):
        raise ValueError("Sample identifiers must be nonempty and unique")
    denied = sorted({accession for accession in scorer.DENY
                     if accession in str(args.counts).upper() or any(accession in str(sample).upper() for sample in samples)})
    if denied:
        raise ValueError(f"Input references release-excluded sources: {denied}")
    result = score_profiles(counts, genes, args.units, args.transport, args.reference)
    args.out.mkdir(parents=True)
    np.savez_compressed(args.out / "profiles.npz", profile=result["profile"], region_key=result["region_key"],
                        sample_id=samples.astype(str), form=FORM, transport=args.transport)
    result["report"]["sample_id"] = samples.tolist()
    result["report"]["input_sha256"] = sha256(args.counts)
    (args.out / "coverage.json").write_text(json.dumps(result["report"], indent=2) + "\n")
    print(json.dumps({"output": str(args.out), "n_samples": len(samples),
                      "transport": args.transport, "coverage": result["report"]["coverage"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, help="Genes-by-samples TSV, first column stable gene IDs")
    parser.add_argument("--units", choices=("counts", "deposited_estimates"))
    parser.add_argument("--transport", choices=("raw", "marginal", "reference_mass"), default="raw")
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--out", type=Path, help="New directory for profiles.npz and coverage.json")
    parser.add_argument("--export-reference", type=Path, help="Create a new fixed training-reference directory")
    main(parser.parse_args())
