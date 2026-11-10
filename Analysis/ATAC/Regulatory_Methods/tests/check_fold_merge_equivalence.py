"""Compute-only full five-fold/shard equivalence and rejection checks."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd
from scipy import sparse


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute allocation required")
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / "scripts"))
    import network_compare as candidate
    import merge_network_folds as merger
    reference = load(root / "results/20261023/network_full_source_v3/scripts/network_compare.py", "reference_fold_fixture")
    args.output.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(42)
    counts = sparse.csr_matrix(rng.poisson(5, size=(30, 12)))
    meta = pd.DataFrame(dict(donor_id=np.repeat([f"d{i}" for i in range(10)], 3), cohort="GSE296875", lineage="hepatocyte"))
    peaks = pd.DataFrame(dict(peak_id=[f"p{i}" for i in range(12)], chrom="chr1", start0=np.arange(12)*2000, end=np.arange(12)*2000+500))
    genes = pd.DataFrame([dict(gene_id=f"g{i}", chrom="chr1", tss0=2000*i, strand="+" if i % 2 else "-") for i in range(1, 13)])
    genes = pd.concat([genes, pd.DataFrame([dict(gene_id="missing", chrom="chrX", tss0=0, strand="-")])], ignore_index=True)
    rna = sparse.csr_matrix(rng.poisson(5, size=(30, 12)))
    rna_genes = pd.DataFrame(dict(gene_id=[f"g{i}" for i in range(1, 13)]))
    graph = pd.DataFrame([dict(peak1=f"p{i}", peak2=f"p{j}", coaccess=(i+j+1)/24) for i in range(12) for j in range(i+1, 12)])
    settings = dict(cicero_rscript="/usr/bin/true", bootstrap_reps=100, max_genes=0, outer_folds=5, inner_folds=3)
    name = "GSE296875_hepatocyte"
    def fake_cicero(*pos, **kwargs):
        stage = pos[2]
        stage.mkdir()
        (stage / "timing.json").write_text(json.dumps(dict(directed_edges=len(graph))))
        (stage / "run.log").write_text('Models with errors:  0\nwindow_models_complete\n')
        (stage / "counts.mtx").write_text("fixture")
        graph.to_csv(stage / "links.tsv", sep="\t", index=False)
        return graph
    for module in [reference, candidate]:
        module.cicero_scores = fake_cicero
    genes.to_csv(args.output / "genes.tsv", sep="\t", index=False)
    def run(module, destination, local_settings):
        destination.mkdir()
        module.fit_group(counts, meta, peaks, genes, destination / name, local_settings, 42, rna, rna_genes)
        (destination / "run_settings.json").write_text(json.dumps(local_settings))
        for filename in [f"{name}_metacell_membership.tsv", f"{name}_eligibility.tsv", "condition_exclusions.tsv"]:
            (destination / filename).write_text("fixture\n")
        snapshot = destination / "source_snapshot"
        snapshot.mkdir()
        for filename in ["network_cicero.R", "network_directnet_native.R"]:
            shutil.copyfile(root / "scripts" / filename, snapshot / filename)
        shutil.copyfile(module.__file__, snapshot / "network_compare.py")
        (destination / "source_sha256.json").write_text(json.dumps({str(args.output / "genes.tsv"): merger.sha(args.output / "genes.tsv")}))
    monolithic = args.output / "monolithic"
    run(reference, monolithic, {**settings, "readout_workers": 1})
    shards = []
    for fold in range(5):
        shard = args.output / f"shard{fold}"
        run(candidate, shard, {**settings, "readout_workers": 6, "selected_outer_folds": [fold]})
        shards.append(shard)
    merged = args.output / "merged"
    merger.merge(shards, [], merged)
    scientific = ["predictions.tsv", "links.tsv", "coverage.tsv", "donor_errors.tsv", "comparisons.tsv", "fold_assignments.tsv", "method_status.tsv", "method_execution_summary.tsv", "failures.tsv"]
    for filename in scientific:
        a, b = monolithic / name / filename, merged / name / filename
        if a.read_bytes() != b.read_bytes():
            raise AssertionError(f"Nonidentical merged table: {filename}")
    hybrid = args.output / "legacy_merged"
    merger.merge(shards[1:], [(monolithic, 0)], hybrid)
    for filename in scientific:
        if (hybrid / name / filename).read_bytes() != (monolithic / name / filename).read_bytes():
            raise AssertionError(f"Nonidentical legacy merged table: {filename}")
    rejected = []
    def reject(label, function):
        try:
            function()
        except (ValueError, FileNotFoundError):
            rejected.append(label)
        else:
            raise AssertionError(f"Invalid merge accepted: {label}")
    reject("missing_fold", lambda: merger.merge(shards[1:], [], args.output / "bad_missing"))
    reject("duplicate_fold", lambda: merger.merge(shards + [shards[0]], [], args.output / "bad_duplicate"))
    path = shards[0] / name / "coverage.tsv"
    saved = path.read_bytes()
    path.write_bytes(saved.splitlines(keepends=True)[0] + b"".join(saved.splitlines(keepends=True)[2:]))
    reject("missing_gene", lambda: merger.merge(shards, [], args.output / "bad_gene"))
    path.write_bytes(saved)
    path = shards[0] / name / "fold_contract.json"
    saved = path.read_bytes()
    contract = json.loads(saved)
    contract["counts_sha256"] = "different"
    path.write_text(json.dumps(contract))
    reject("different_counts", lambda: merger.merge(shards, [], args.output / "bad_counts"))
    path.write_bytes(saved)
    path = shards[0] / name / "cicero_fold0_inner0/run.log"
    saved = path.read_bytes()
    path.write_text("Models with errors:  1\nwindow_models_complete\n")
    reject("cicero_failure", lambda: merger.merge(shards, [], args.output / "bad_cicero"))
    path.write_bytes(saved)
    for invalid in [[], [0, 0], [-1], [5], [True], "0"]:
        reject(f"invalid_selector_{invalid}", lambda: candidate.selected_folds(dict(selected_outer_folds=invalid), range(5)))
    (args.output / "verification.json").write_text(json.dumps(dict(status="passed", exact_scientific_tables=scientific,
        monolithic_workers=1, shard_workers=6, folds=5, legacy_snapshot_fold=0, rejected=rejected,
        candidate_source_sha256=merger.sha(Path(candidate.__file__))), indent=2) + "\n")


if __name__ == "__main__":
    main()
