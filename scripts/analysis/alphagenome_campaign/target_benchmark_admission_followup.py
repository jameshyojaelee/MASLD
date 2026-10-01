#!/usr/bin/env python3
"""Fixed distance-only feasibility, activity-source inventory, independent F5 counts.

No model fitting, saturation reads, threshold search, allele inference or liver
validation. F5 is independently reconstructed using only original source tables.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import urllib.request

import numpy as np
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "GWAS/finemapping/results/alphagenome_campaign/week1-20260915/data"
CTX = ROOT / "Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1"
PROGRAMS = ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv"
UNIVERSE = ROOT / "GWAS/finemapping/results/alphagenome_atlas/p4-saturation-20260909T191917Z/tables/saturation_universe.tsv.gz"


def rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def dump(path, result):
    path.write_text(json.dumps(result, indent=2) + "\n")


def f5_independent(out, producer):
    # Independent row-wise implementation; no imports from the producer.
    permitted = set()
    for r in rows(CTX / "consensus_peak_manifest.tsv"):
        if r["blacklist_overlap"].upper() != "TRUE":
            permitted.add((r["lineage"], f"{r['chrom']}:{r['start0']}-{r['end']}"))
    source_relations = set()
    for r in rows(CTX / "da/da_peak_results.tsv.gz"):
        if (r["lineage"], r["peak_coordinate"]) not in permitted:
            continue
        for gene in r["promoter_genes"].replace(";", ",").split(","):
            if gene not in {"", "NA", "nan", "none", "None"}:
                source_relations.add((gene, r["peak_coordinate"]))
    registered = {r["region_key"] for r in rows(UNIVERSE) if r["universe"] == "U2_snatac"}
    assert len(registered) == 6930
    membership = collections.defaultdict(set)
    for r in rows(PROGRAMS):
        membership[r["program_uid"]].add(r["canonical_gene"])
    assert len(membership) == 117
    by_gene = collections.defaultdict(set)
    for gene, region in source_relations:
        by_gene[gene].add(region)
    reconstructed_background = {}
    for gene, regions in by_gene.items():
        reconstructed_background[gene] = (len(regions), sum(r in registered for r in regions), sum(r not in registered for r in regions))
    reconstructed_programs = {}
    for uid, genes in membership.items():
        regions = {region for gene in genes for region in by_gene.get(gene, ())}
        n_registered = sum(region in registered for region in regions)
        reconstructed_programs[uid] = (len(genes), len(regions), n_registered, len(regions) - n_registered, n_registered >= 200)
    reported_background = {r["canonical_gene"]: tuple(int(r[k]) for k in ("eligible_regions", "registered_regions", "missing_regions"))
                           for r in rows(producer / "original_background_coverage.tsv")}
    reported_programs = {r["program_uid"]: tuple(int(r[k]) for k in ("program_genes", "eligible_regions", "registered_regions", "missing_regions"))
                        + (r["original_minimum200_possible_on_registered_universe"].lower() == "true",)
                        for r in rows(producer / "original_program_coverage.tsv")}
    missing = {region for _, region in source_relations if region not in registered}
    reported_missing = {r["region_key"] for r in rows(producer / "missing_background_regions.tsv.gz")}
    comparisons = {"all_gene_rows_identical": reconstructed_background == reported_background,
                   "all_program_rows_identical": reconstructed_programs == reported_programs,
                   "all_missing_region_identities_identical": missing == reported_missing}
    result = {"comparisons": comparisons, "background_genes": len(by_gene),
              "background_unique_regions": len({r for _, r in source_relations}),
              "missing_unique_regions": len(missing),
              "background_genes_with_any_missing_region": sum(v[2] > 0 for v in reconstructed_background.values()),
              "programs": len(membership), "programs_with_200_registered_regions": sum(v[4] for v in reconstructed_programs.values()),
              "registered_program_region_range": [min(v[2] for v in reconstructed_programs.values()), max(v[2] for v in reconstructed_programs.values())],
              "no_saturation_values_read": True, "registered_is_not_retrieved_or_intact": True,
              "F5_complete": False, "producer_output": str(producer),
              "source_program_sha256": sha(PROGRAMS), "source_universe_sha256": sha(UNIVERSE)}
    dump(out / "F5_independent_reconstruction.json", result)
    if not all(comparisons.values()):
        raise ValueError("Independent original F5 metadata reconstruction disagrees")


def average_precision(y, score, weights=None):
    if weights is None:
        weights = np.ones(len(y))
    if weights.sum() == 0 or np.dot(y, weights) == 0:
        return None
    return float(average_precision_score(y, score, sample_weight=weights))


def metric(part, label):
    y = np.array([int(r[label]) for r in part])
    score = -np.log10(np.array([abs(float(r["distance_to_TSS_source"])) for r in part]))
    chromosomes = np.array([r["element_chrom_source"] for r in part])
    levels = np.unique(chromosomes)
    observed = average_precision(y, score)
    rng = np.random.default_rng(20260916)
    boot = []
    for _ in range(500):
        weights_by_chr = collections.Counter(rng.choice(levels, len(levels), replace=True))
        weights = np.array([weights_by_chr[c] for c in chromosomes])
        value = average_precision(y, score, weights)
        if value is not None:
            boot.append(value)
    groups = collections.defaultdict(list)
    for r, label_value, s in zip(part, y, score):
        groups[(r["source_asset"], r["cell_type"], r["element_chrom_source"],
                r["element_start_source"], r["element_end_source"])].append((r["gene_id_source"], int(label_value), s))
    multitarget = [v for v in groups.values() if len({g for g, _, _ in v}) >= 2]
    top_precision = []
    detected = 0
    for group in multitarget:
        if len({g for g, _, _ in group}) != len(group):
            raise ValueError("Repeated same-source target measurement requires explicit collapse before baseline")
        maximum = max(s for _, _, s in group)
        top = [label_value for _, label_value, s in group if s == maximum]
        top_precision.append(float(np.mean(top)))
        detected += int(any(label_value for _, label_value, _ in group))
    return {"endpoint": label, "pairs": len(part), "positive_pairs": int(y.sum()),
            "chromosomes": len(levels), "source_elements": len(groups),
            "constant_score_average_precision": float(y.mean()), "distance_average_precision": observed,
            "distance_average_precision_chr_bootstrap_95CI": [float(x) for x in np.quantile(boot, [.025, .975])] if boot else None,
            "bootstrap_replicates": 500, "bootstrap_finite_replicates": len(boot),
            "multitarget_elements": len(multitarget), "multitarget_elements_with_detected_target": detected,
            "top_ranked_measured_target_precision_tie_average": float(np.mean(top_precision)) if top_precision else None,
            "multitarget_population": "defined by >=2 eligible measured genes irrespective of outcomes; elements lacking detected targets retained",
            "interpretation": "fixed_unfitted_distance_baseline;source_native_distance;not_liver_validation_or_external_model_generalization"}


def target_feasibility(admission, out):
    grouped = collections.defaultdict(list)
    for asset in ("klann", "morrisv1", "morrisv2", "xie", "dc_tap_source"):
        for r in rows(admission / f"{asset}.normalized.tsv.gz"):
            if r["binary_endpoint_eligible"] != "True":
                continue
            if int(r["same_tuple_source_row_multiplicity"]) != 1:
                raise ValueError("Duplicate endpoint rows require explicit handling")
            grouped[(asset, r["cell_type"], "all_geometry_uniform_power")].append(r)
            if asset == "dc_tap_source" and r["source_random_distal_design"].upper() == "TRUE":
                grouped[(asset, r["cell_type"], "native_random_distal_design")].append(r)
    outputs = []
    for (asset, cell, population), part in sorted(grouped.items()):
        for label in ("expression_decrease_binary_label", "any_expression_change_binary_label"):
            outputs.append(dict(source_asset=asset, cell_type=cell, population=population, **metric(part, label)))
    # Scientific invariants: ties do not depend on row order; a constant score
    # has AP equal to observed prevalence; no-label-positive population is unknown.
    assert average_precision(np.array([1, 0, 1, 0]), np.ones(4)) == .5
    assert average_precision(np.array([0, 1, 0, 1]), np.ones(4)) == .5
    assert average_precision(np.zeros(4), np.arange(4)) is None
    dump(out / "target_distance_feasibility.json", {"results": outputs,
         "AP_definition": "sklearn noninterpolated average precision;ties evaluated at shared thresholds",
         "uncertainty": "500 chromosome-resampled descriptive bootstraps;no significance tests or multiplicity claim",
         "assay_units": "binary endpoints kept separate;no pooled expression-effect regression",
         "activity_comparison": "not_run_pending_all_candidate_source_identity_and_coverage",
         "invariants": {"tied_score_prevalence": True, "zero_positive_undefined_AP": True}})


def activity_inventory(out, manifest):
    asset = json.loads(manifest.read_text())
    dest = out / asset["filename"]
    expected = asset["bytes"]
    git_hash = hashlib.sha1(f"blob {expected}\0".encode())
    count = 0
    with urllib.request.urlopen(asset["url"], timeout=60) as source, dest.open("xb") as handle:
        for block in iter(lambda: source.read(1 << 20), b""):
            count += len(block)
            if count > expected:
                raise ValueError("activity source exceeded fixed size")
            git_hash.update(block)
            handle.write(block)
    verified = count == expected and git_hash.hexdigest() == asset["git_blob_sha1"]
    receipt = dict(asset, observed_bytes=count, observed_git_blob_sha1=git_hash.hexdigest(), sha256=sha(dest), verified=verified)
    dump(out / "activity_source_receipt.json", receipt)
    if not verified:
        raise ValueError("Activity source blob verification failed")
    n = 0
    columns = []
    missing = collections.Counter()
    examples = []
    for row in rows(dest):
        n += 1
        if not columns:
            columns = list(row)
        for k, v in row.items():
            if v in {"", "NA", "NaN", "nan"}:
                missing[k] += 1
        if n <= 2:
            examples.append({k: v for k, v in row.items() if any(t in k.lower() for t in ["chrom", "start", "end", "gene", "cell", "name", "dhs", "dnase", "rpm"])})
    dump(out / "activity_source_inventory.json", {"rows": n, "columns": columns, "missing_by_column": dict(missing),
         "identity_activity_examples_only": examples, "status": "inventory_only_no_cross_version_or_resized_region_join_assumed",
         "next_required": "verify original-to-resized region mapping and cell-type identity on all source candidates before any activity baseline"})


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--admission", type=Path, required=True)
    p.add_argument("--f5-producer", type=Path, required=True)
    p.add_argument("--activity-manifest", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    args.out.mkdir(parents=True, exist_ok=False)
    target_feasibility(args.admission, args.out)
    f5_independent(args.out, args.f5_producer)
    activity_inventory(args.out, args.activity_manifest)
    dump(args.out / "receipt.json", {"status": "complete", "job_id": os.environ["SLURM_JOB_ID"],
         "script_sha256": sha(Path(__file__)), "python": sys.version, "no_model_fit": True,
         "no_saturation_inference_or_raw_reads": True, "primary_liver_target_validation": False})
    print("distance baseline, independent F5 inventory and activity-source inspection complete", flush=True)


if __name__ == "__main__":
    main()
