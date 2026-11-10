#!/usr/bin/env python3
"""Merge disjoint finished donor-held-out folds; recompute donor inference once.

Input paths are network CLI analysis roots, each containing exactly one group.
Legacy snapshots require the explicitly reviewed frozen-v3 numerical source.
Never opens an existing output for writing.
"""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
import re
import shutil

import pandas as pd

import network_compare as core


LEGACY_SOURCE = "f0d9de80dca12c7a97ef5dad3cf879f18548fc522e63ebd577142a14655bd26d"
TABLES = ["predictions.tsv", "links.tsv", "coverage.tsv", "gene_timing.tsv"]


def read_json(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def group_path(root):
    groups = sorted(p.parent for p in root.glob("*/fold_assignments.tsv"))
    if len(groups) != 1:
        raise ValueError(f"Expected exactly one group in {root}")
    return groups[0]


def selected_rows(path, fold):
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if int(row["outer_fold"]) == fold:
                yield row


def verify_cicero(group, fold):
    for tag in ["outer", "inner0", "inner1", "inner2"]:
        stage = group / f"cicero_fold{fold}_{tag}"
        timing = read_json(stage / "timing.json")
        log = (stage / "run.log").read_text()
        errors = re.findall(r"Models with errors:\s*(\d+)", log)
        if errors != ["0"] or "window_models_complete" not in log or timing.get("directed_edges", 0) <= 0:
            raise ValueError(f"Cicero stage not successfully completed: {stage}")
        if not (stage / "links.tsv").is_file():
            raise ValueError(f"Missing Cicero graph: {stage}")


def verify_common(root, reference):
    for authority in [root, group_path(root), root.parent, root.parent.parent]:
        if any(authority.glob("INVALIDATED*.json")):
            raise ValueError(f"Invalidated source: {authority}")
    current = group_path(root)
    expected = group_path(reference)
    if current.name != expected.name:
        raise ValueError("Different cohort-lineage groups")
    a, b = read_json(root / "run_settings.json"), read_json(reference / "run_settings.json")
    if core.scientific_settings(a) != core.scientific_settings(b):
        raise ValueError("Different scientific settings")
    for relative in [f"{current.name}_metacell_membership.tsv", f"{current.name}_eligibility.tsv", "condition_exclusions.tsv"]:
        if sha(root / relative) != sha(reference / relative):
            raise ValueError(f"Different input roster: {relative}")
    for relative in ["fold_assignments.tsv", "method_status.tsv", "dependencies.json"]:
        if sha(current / relative) != sha(expected / relative):
            raise ValueError(f"Different fold/method inputs: {relative}")
    # Python changed only through the explicitly reviewed legacy exception below.
    for name in ["network_cicero.R", "network_directnet_native.R"]:
        if sha(root / "source_snapshot" / name) != sha(reference / "source_snapshot" / name):
            raise ValueError(f"Different method source: {name}")
    source_a = read_json(root / "source_sha256.json")
    source_b = read_json(reference / "source_sha256.json")
    annotation_a = {Path(k).name: v for k, v in source_a.items() if Path(k).suffix not in {".py", ".R", ".json"}}
    annotation_b = {Path(k).name: v for k, v in source_b.items() if Path(k).suffix not in {".py", ".R", ".json"}}
    if annotation_a != annotation_b or not annotation_a:
        raise ValueError("Different or missing annotation source identity")


def qualify_snapshot(root, fold, reference, receipt):
    """Qualify immutable v3 fold0 rows while subsequent folds keep running.

    The reference need only have finished input preparation and its contract.
    Hashes describe selected-fold rows, since later-fold streams may grow.
    """
    if fold != 0:
        raise ValueError("Early legacy qualification is restricted to original fold0")
    verify_common(reference, reference)
    verify_common(root, reference)
    group, ref_group = group_path(root), group_path(reference)
    contract = read_json(ref_group / "fold_contract.json")
    if sha(reference / "source_snapshot" / "network_compare.py") != contract["source_sha256"]:
        raise ValueError("Reference source differs from its input contract")
    if sha(ref_group / "fold_assignments.tsv") != contract["fold_assignments_sha256"]:
        raise ValueError("Reference donor folds differ from their contract")
    if sha(root / "source_snapshot" / "network_compare.py") != LEGACY_SOURCE:
        raise ValueError("Unreviewed legacy numerical source")
    if not (group / "cicero_fold1_outer" / "counts.mtx").exists():
        raise ValueError("Legacy fold has not crossed its stream-flush boundary")
    verify_cicero(group, fold)
    coverage = list(selected_rows(group / "coverage.tsv", fold))
    if [r["gene_id"] for r in coverage] != contract["gene_ids"]:
        raise ValueError("Missing, duplicated or reordered gene coverage")
    timings = list(selected_rows(group / "gene_timing.tsv", fold))
    if [r["gene_id"] for r in timings] != [r["gene_id"] for r in coverage if r["status"] == "testable"]:
        raise ValueError("Missing, duplicated or unflushed completed-gene timing records")
    timing_by_gene = {r["gene_id"]: r for r in timings}
    positions = {gene: i for i, gene in enumerate(contract["gene_ids"])}
    metadata = pd.read_csv(group / "fold_assignments.tsv", sep="\t", dtype=str, keep_default_na=False)
    held_out = metadata[metadata.outer_fold == str(fold)].groupby("donor_id").size().to_dict()
    method_status = pd.read_csv(group / "method_status.tsv", sep="\t", keep_default_na=False)
    methods = set(method_status[method_status.status == "available"].method)
    methods.update(["intercept", "promoter_only"])
    counts, hashes = {}, {}
    for filename, count_column in [("predictions.tsv", "n_prediction_rows"), ("links.tsv", "n_link_rows")]:
        digest = hashlib.sha256()
        total, previous, observed = 0, -1, set()
        rows = selected_rows(group / filename, fold)
        for gene, block in itertools.groupby(rows, key=lambda row: row["gene_id"]):
            position = positions.get(gene, -1)
            if position <= previous or gene not in timing_by_gene:
                raise ValueError(f"Unknown, duplicated or reordered gene in {filename}")
            previous = position
            observed.add(gene)
            keys, peak_keys, gene_count = set(), set(), 0
            for row in block:
                if row["method"] not in methods:
                    raise ValueError("Unknown prediction/link method")
                if filename == "predictions.tsv":
                    key = (row["endpoint"], row["method"], row["donor_id"])
                    if row["endpoint"] not in contract["planned_endpoints"] or int(row["n_metacells"]) != held_out.get(row["donor_id"]):
                        raise ValueError("Prediction donor or endpoint lies outside its held-out fold")
                else:
                    key = (row["method"], int(row["rank"]))
                    peak_key = (row["method"], row["peak_id"])
                    if peak_key in peak_keys or int(row["rank"]) < 1:
                        raise ValueError("Duplicate peak or invalid link rank")
                    peak_keys.add(peak_key)
                if key in keys:
                    raise ValueError(f"Duplicate scientific row in {filename}")
                keys.add(key)
                digest.update((json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode())
                gene_count += 1
            if gene_count != int(timing_by_gene[gene][count_column]):
                raise ValueError(f"Missing or extra gene rows in {filename}")
            total += gene_count
        expected = {gene for gene, row in timing_by_gene.items() if int(row[count_column]) > 0}
        if observed != expected:
            raise ValueError(f"Missing completed-gene rows in {filename}")
        counts[filename], hashes[filename] = total, digest.hexdigest()
    for filename in ["coverage.tsv", "gene_timing.tsv"]:
        digest = hashlib.sha256()
        total = 0
        for row in selected_rows(group / filename, fold):
            digest.update((json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode())
            total += 1
        counts[filename], hashes[filename] = total, digest.hexdigest()
    result = dict(status="qualified_snapshot", outer_fold=fold, source_root=str(root.resolve()),
        reference_root=str(reference.resolve()), legacy_source_sha256=LEGACY_SOURCE,
        reference_contract_sha256=sha(ref_group / "fold_contract.json"),
        inference_complete=False, scientific_adoption=False, n_genes=len(coverage),
        row_counts=counts, selected_fold_row_sha256=hashes,
        qualification="reviewed_v3_source_equal_settings_annotation_rosters_folds_and_finished_fold_rows")
    receipt.parent.mkdir(parents=True, exist_ok=True)
    with receipt.open("x") as handle:
        handle.write(json.dumps(result, indent=2) + "\n")
    return result


def merge(shards, snapshots, output):
    if not shards:
        raise ValueError("At least one contracted completed shard is required")
    reference = shards[0]
    ref_group = group_path(reference)
    contract = read_json(ref_group / "fold_contract.json")
    available = contract["available_folds"]
    if available != list(range(5)) or int(contract["settings"].get("inner_folds", 3)) != 3:
        raise ValueError("Merger requires the prespecified five outer/three inner folds")
    sources = {}
    for root in shards:
        verify_common(root, reference)
        group = group_path(root)
        if read_json(group / "fold_contract.json") != contract:
            raise ValueError("Different aggregate inputs, source or analysis contract")
        if sha(root / "source_snapshot" / "network_compare.py") != contract["source_sha256"]:
            raise ValueError("Shard numerical source differs from its contract")
        if sha(group / "fold_assignments.tsv") != contract["fold_assignments_sha256"]:
            raise ValueError("Shard donor folds differ from their contract")
        marker = read_json(group / "shard_completion.json")
        if marker.get("status") != "completed" or marker.get("all_outer_folds") is not False:
            raise ValueError("Unfinished or incorrectly labeled shard")
        coverage_folds = set()
        with (group / "coverage.tsv").open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                coverage_folds.add(int(row["outer_fold"]))
        if coverage_folds != set(marker["selected_outer_folds"]):
            raise ValueError("Shard coverage folds differ from its completion marker")
        for fold in marker["selected_outer_folds"]:
            finished = read_json(group / f"fold{fold}_completed.json")
            if finished.get("status") != "completed_fold" or finished.get("failures") or finished.get("n_genes") != len(contract["gene_ids"]):
                raise ValueError("Incomplete or failed fold")
            if fold in sources:
                raise ValueError("Duplicate fold")
            sources[fold] = (root, group, False)
    for root, fold in snapshots:
        verify_common(root, reference)
        group = group_path(root)
        if sha(root / "source_snapshot" / "network_compare.py") != LEGACY_SOURCE:
            raise ValueError("Unreviewed legacy numerical source")
        # v3 flushes all streams before beginning the next fold's outer stage.
        if not (group / f"cicero_fold{fold + 1}_outer" / "counts.mtx").exists():
            raise ValueError("Legacy fold has not crossed its stream-flush boundary")
        if fold in sources:
            raise ValueError("Duplicate fold")
        sources[fold] = (root, group, True)
    if sorted(sources) != available:
        raise ValueError("Missing or unexpected outer folds")
    for fold, (_, group, _) in sources.items():
        verify_cicero(group, fold)
        coverage = list(selected_rows(group / "coverage.tsv", fold))
        observed = [r["gene_id"] for r in coverage]
        if observed != contract["gene_ids"]:
            raise ValueError("Missing, duplicated or reordered gene coverage")
        timings = list(selected_rows(group / "gene_timing.tsv", fold))
        if [r["gene_id"] for r in timings] != [r["gene_id"] for r in coverage if r["status"] == "testable"]:
            raise ValueError("Missing, duplicated or unflushed completed-gene timing records")
    output.mkdir(parents=True, exist_ok=False)
    target = output / ref_group.name
    target.mkdir()
    for relative in ["run_settings.json", "source_sha256.json", "condition_exclusions.tsv",
                     f"{ref_group.name}_metacell_membership.tsv", f"{ref_group.name}_eligibility.tsv"]:
        shutil.copyfile(reference / relative, output / relative)
    run_settings = read_json(output / "run_settings.json")
    run_settings.pop("selected_outer_folds", None)
    (output / "run_settings.json").write_text(json.dumps(run_settings, indent=2) + "\n")
    shutil.copytree(reference / "source_snapshot", output / "source_snapshot")
    for name in ["fold_assignments.tsv", "method_status.tsv", "dependencies.json", "memory_projection.json", "fold_contract.json"]:
        shutil.copyfile(ref_group / name, target / name)
    accumulator = core.ComparisonAccumulator(contract["planned_endpoints"], ["distance", "directnet_adapted", "cicero"])
    if contract["settings"].get("external_stage_links"):
        raise ValueError("External-stage merger is outside this primary comparison")
    counts, scored_folds = {}, {}
    fold_metadata = pd.read_csv(ref_group / "fold_assignments.tsv", sep="\t", dtype=str, keep_default_na=False)
    allowed_donors = {fold: fold_metadata[fold_metadata.outer_fold == str(fold)].groupby("donor_id").size().to_dict() for fold in available}
    with (target / "predictions.tsv").open("x", newline="") as handle:
        columns = next(csv.reader((ref_group / "predictions.tsv").open(), delimiter="\t"))
        writer = csv.DictWriter(handle, columns, delimiter="\t")
        writer.writeheader()
        for fold in available:
            group = sources[fold][1]
            expected_rows = {r["gene_id"]: int(r["n_prediction_rows"]) for r in selected_rows(group / "gene_timing.tsv", fold) if int(r["n_prediction_rows"]) > 0}
            observed_genes = set()
            last_gene = -1
            gene_positions = {gene: i for i, gene in enumerate(contract["gene_ids"])}
            rows = selected_rows(group / "predictions.tsv", fold)
            for gene, block in itertools.groupby(rows, key=lambda r: r["gene_id"]):
                position = gene_positions.get(gene, -1)
                if position <= last_gene:
                    raise ValueError("Unknown, duplicated or reordered prediction gene")
                last_gene = position
                block = list(block)
                if len(block) != expected_rows.get(gene):
                    raise ValueError("Missing or extra gene prediction rows")
                observed_genes.add(gene)
                keys = [(r["endpoint"], r["method"], r["donor_id"]) for r in block]
                if len(keys) != len(set(keys)):
                    raise ValueError("Duplicate donor prediction")
                numeric = []
                for row in block:
                    if int(row["n_metacells"]) != allowed_donors[fold].get(row["donor_id"]):
                        raise ValueError("Prediction donor is outside its held-out fold or has wrong coverage")
                    item = dict(row)
                    item["standardized_squared_error"] = float(item["standardized_squared_error"])
                    item["n_metacells"] = int(item["n_metacells"])
                    numeric.append(item)
                    method = row["method"]
                    counts[method] = counts.get(method, 0) + 1
                    scored_folds.setdefault(method, set()).add(fold)
                accumulator.append(numeric)
                writer.writerows(block)
            if observed_genes != set(expected_rows):
                raise ValueError("Missing completed-gene predictions")
    for name in TABLES[1:]:
        with (target / name).open("x", newline="") as handle:
            columns = next(csv.reader((ref_group / name).open(), delimiter="\t"))
            writer = csv.DictWriter(handle, columns, delimiter="\t")
            writer.writeheader()
            for fold in available:
                rows = selected_rows(sources[fold][1] / name, fold)
                expected_links = {r["gene_id"]: int(r["n_link_rows"]) for r in selected_rows(sources[fold][1] / "gene_timing.tsv", fold)} if name == "links.tsv" else None
                actual_links = {}
                previous = -1
                seen = set()
                seen_peaks = set()
                for row in rows:
                    position = gene_positions.get(row["gene_id"], -1)
                    key = (row["gene_id"], row.get("method"), row.get("rank"))
                    peak_key = (row["gene_id"], row.get("method"), row.get("peak_id"))
                    if position < 0 or position < previous or key in seen or (name == "links.tsv" and peak_key in seen_peaks):
                        raise ValueError("Reordered or duplicate fold table rows")
                    if position != previous:
                        seen.clear()
                        seen_peaks.clear()
                    seen.add(key)
                    seen_peaks.add(peak_key)
                    previous = position
                    if expected_links is not None:
                        actual_links[row["gene_id"]] = actual_links.get(row["gene_id"], 0) + 1
                    writer.writerow(row)
                if expected_links is not None and any(actual_links.get(gene, 0) != total for gene, total in expected_links.items()):
                    raise ValueError("Missing or extra completed-gene links")
    statuses = pd.read_csv(target / "method_status.tsv", sep="\t", keep_default_na=False)
    summary = [dict(method=r.method, eligibility=r.status, n_prediction_rows=counts.get(r.method, 0),
                    n_scored_outer_folds=len(scored_folds.get(r.method, set())), n_failed_stages=0,
                    execution_status="completed" if counts.get(r.method, 0) else "not_scored", reason=r.reason)
               for r in statuses.itertuples()]
    pd.DataFrame(summary).to_csv(target / "method_execution_summary.tsv", sep="\t", index=False)
    pd.DataFrame(columns=["outer_fold", "method", "stage", "reason"]).to_csv(target / "failures.tsv", sep="\t", index=False)
    comparisons, donor_errors = accumulator.finalize(contract["seed"], int(contract["settings"].get("bootstrap_reps", 2000)))
    comparisons.to_csv(target / "comparisons.tsv", sep="\t", index=False)
    donor_errors.to_csv(target / "donor_errors.tsv", sep="\t", index=False)
    metadata = pd.read_csv(target / "fold_assignments.tsv", sep="\t", keep_default_na=False)
    comparisons["cohort"], comparisons["lineage"] = metadata.cohort.iloc[0], metadata.lineage.iloc[0]
    comparisons["lineage_state"] = ";".join(sorted(metadata.lineage_state.unique())) if "lineage_state" in metadata else "source_label"
    comparisons["source_lineage"] = ";".join(sorted(metadata.source_lineage.unique())) if "source_lineage" in metadata else metadata.lineage.iloc[0]
    comparisons["holm_pvalue"] = core.holm(comparisons.pvalue.to_numpy())
    comparisons["provisional_predictive_improvement"] = (comparisons.ci_low > 0) & (comparisons.holm_pvalue < 0.05)
    comparisons["predictive_improvement"], comparisons["family_complete"] = False, False
    comparisons["family"] = "provisional_within_invocation_all_planned_contrasts_global_pooling_required"
    comparisons.to_csv(output / "comparisons.tsv", sep="\t", index=False)
    expected = {"distance", "directnet_adapted", "cicero"}
    qualified = all(len(scored_folds.get(method, set())) == 5 for method in expected)
    (output / "completion.json").write_text(json.dumps(dict(status="completed" if qualified else "partial",
        method_fold_qualification=[dict(group=ref_group.name, all_methods_all_outer_folds=qualified)],
        interpretation="candidate_accessibility_context_and_internal_prediction", pilot=False,
        rna_available=contract["rna_sha256"] is not None, scientific_adoption=False), indent=2) + "\n")
    receipt = dict(status="completed" if qualified else "partial", folds=available, sources={str(f): dict(root=str(v[0]), legacy_snapshot=v[2]) for f, v in sources.items()},
                   legacy_input_identity="reviewed_v3_source_equal_settings_annotation_rosters_and_folds" if snapshots else None)
    (output / "merge_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard", type=Path, action="append", required=True)
    parser.add_argument("--snapshot", action="append", default=[], help="ANALYSIS_ROOT:FOLD")
    parser.add_argument("--qualify-snapshot", action="store_true", help="Qualify one finished v3 fold against one prepared shard; --output is a new JSON receipt")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    core.require_compute()
    snapshots = [(Path(value.rsplit(":", 1)[0]), int(value.rsplit(":", 1)[1])) for value in args.snapshot]
    if args.qualify_snapshot:
        if len(args.shard) != 1 or len(snapshots) != 1:
            raise ValueError("Snapshot qualification needs exactly one prepared shard and one snapshot")
        qualify_snapshot(*snapshots[0], args.shard[0], args.output)
        return
    merge(args.shard, snapshots, args.output)


if __name__ == "__main__":
    main()
