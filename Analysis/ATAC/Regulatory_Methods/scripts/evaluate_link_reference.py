#!/usr/bin/env python3
"""Descriptive recovery of supplied qualified references with matched unlabeled edges.

References must contain reference_type/source_id/target_id/cohort/lineage,
reference_id/source_uri/qualification_status. Types are cis_link, TF_binding,
and perturbation, and qualification_status must be 'qualified'. A supplied,
outcome-blind candidate universe needs the first five identity columns plus
distance_bin/observability_bin/observable. Native network links adapt only to
cis_link; TF evaluation requires explicit TF identities, never motif-name guesses.

Optional priors require the edge identity columns plus method/prior_id, with
reference_id identifying any reused reference source. Prior-exposed positives
are audited and excluded. No prior file means independence is unassessed.
Missing reference edges are unlabeled, never verified negatives; no AUROC/AUPRC.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import socket

import numpy as np
import pandas as pd


TYPES = {"cis_link", "TF_binding", "perturbation"}
IDENTITY = ["reference_type", "cohort", "lineage", "source_id", "target_id"]
GROUP = ["reference_type", "cohort", "lineage", "method", "outer_fold"]


def table(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def require_columns(frame, columns, label):
    if not set(columns).issubset(frame):
        raise ValueError(f"{label} requires {columns}")
    if frame[columns].eq("").any().any():
        raise ValueError(f"{label} has blank identifying/required values")


def validate_reference(frame):
    require_columns(frame, IDENTITY + ["reference_id", "source_uri", "qualification_status"], "Reference")
    if not set(frame.reference_type).issubset(TYPES):
        raise ValueError("Unknown reference type")
    return frame[frame.qualification_status.eq("qualified")].copy()


def adapt_links(path):
    frame = table(path)
    if frame.empty:
        return pd.DataFrame(columns=GROUP + ["source_id", "target_id", "score"])
    # A native accessibility link cannot be converted into a unique TF edge.
    if not set(IDENTITY).issubset(frame):
        require_columns(frame, ["peak_id", "gene_id", "method", "outer_fold", "score"], "Native link")
        if "cohort" not in frame or "lineage" not in frame:
            name = path.parent.name.split("_", 1)
            if len(name) != 2 or not name[0].startswith("GSE"):
                raise ValueError("Native links require cohort/lineage columns or GSE*_lineage parent")
            frame["cohort"], frame["lineage"] = name
        frame["reference_type"] = "cis_link"
        frame["source_id"], frame["target_id"] = frame.peak_id, frame.gene_id
    require_columns(frame, GROUP + ["source_id", "target_id", "score"], "Prediction")
    if not set(frame.reference_type).issubset(TYPES):
        raise ValueError("Unknown prediction reference type")
    if "motif_not_unique_tf" in frame:
        ambiguous = frame.motif_not_unique_tf.str.lower().isin(["true", "1"])
        if (ambiguous & ~frame.reference_type.eq("cis_link")).any():
            raise ValueError("Ambiguous motif compatibility cannot establish a unique TF prediction")
    frame["score"] = pd.to_numeric(frame.score, errors="raise")
    if not np.isfinite(frame.score).all():
        raise ValueError("Prediction scores must be finite")
    if frame.duplicated(GROUP + ["source_id", "target_id"]).any():
        raise ValueError("Duplicate prediction edge within context/method/fold")
    return frame


def top_edges(frame, top_k):
    # Per-target ranking avoids treating gene-specific XGBoost gains as a global scale.
    ordered = frame.sort_values(["target_id", "score", "source_id"], ascending=[True, False, True], kind="stable")
    return set(map(tuple, ordered.groupby("target_id", sort=False).head(top_k)[["source_id", "target_id"]].to_numpy()))


def evaluate_group(predictions, references, universe, priors, rng, top_k, repeats):
    context = {key: str(predictions.iloc[0][key]) for key in GROUP}
    subset = lambda f: f[(f.reference_type == context["reference_type"]) & (f.cohort == context["cohort"]) & (f.lineage == context["lineage"])]
    reference = subset(references).copy()
    candidate = subset(universe).copy()
    result = {**context, "n_qualified_reference_edges": len(reference), "top_k_per_target": top_k,
              "background_repeats": repeats, "status": "not_run", "reason": "no_qualified_reference_for_context"}
    if reference.empty:
        return result, []
    prior = None if priors is None else subset(priors)[lambda x: x.method.eq(context["method"])]
    prior_edges = set() if prior is None else set(map(tuple, prior[["source_id", "target_id"]].to_numpy()))
    prior_sources = set() if prior is None or "reference_id" not in prior else set(prior.reference_id) - {""}
    audit = []
    retained = []
    for row in reference.to_dict("records"):
        state = "unassessed" if prior is None else "exact_edge_overlap" if (row["source_id"], row["target_id"]) in prior_edges else "reference_source_overlap" if row["reference_id"] in prior_sources else "no_recorded_overlap"
        audit.append({**context, **{key: row[key] for key in ["source_id", "target_id", "reference_id", "source_uri"]}, "prior_overlap_state": state})
        if state in {"unassessed", "no_recorded_overlap"}:
            retained.append(row)
    result["n_prior_exposed_reference_edges"] = len(reference) - len(retained)
    result["independence_state"] = "prior_overlap_unassessed" if prior is None else "supplied_prior_manifest_audited"
    if not retained:
        result["reason"] = "all_reference_edges_exposed_to_prior"
        return result, audit
    reference = pd.DataFrame(retained)
    all_positive = set(map(tuple, subset(references)[["source_id", "target_id"]].to_numpy()))
    observable = candidate.observable.str.lower().isin(["true", "1"])
    candidate = candidate[observable]
    joined = reference[IDENTITY + ["reference_id", "source_uri"]].merge(candidate, on=IDENTITY, how="inner", validate="one_to_one")
    result["n_observable_unexposed_reference_edges"] = len(joined)
    if joined.empty:
        result["reason"] = "no_observable_reference_edges"
        return result, audit
    measured_edges = set(map(tuple, candidate[["source_id", "target_id"]].to_numpy()))
    measured_predictions = predictions[[tuple(edge) in measured_edges for edge in predictions[["source_id", "target_id"]].to_numpy()]]
    result["n_prediction_edges_outside_observable_universe"] = len(predictions) - len(measured_predictions)
    selected = top_edges(measured_predictions, top_k)
    # Hold target gene fixed for cis links, and TF fixed for TF-binding/perturbation.
    anchor = "target_id" if context["reference_type"] == "cis_link" else "source_id"
    strata = [anchor, "distance_bin", "observability_bin"]
    background = candidate[["source_id", "target_id"] + [c for c in strata if c not in ["source_id", "target_id"]]].copy()
    background = background[[tuple(edge) not in all_positive for edge in background[["source_id", "target_id"]].to_numpy()]]
    pools = {key: list(map(tuple, value[["source_id", "target_id"]].to_numpy())) for key, value in background.groupby(strata, sort=True)}
    counts, observed, matched, excluded = [], 0, 0, 0
    for key, positive in joined.groupby(strata, sort=True):
        pool = pools.get(key, [])
        if len(pool) < len(positive):
            excluded += len(positive)
            continue
        observed += sum(tuple(edge) in selected for edge in positive[["source_id", "target_id"]].to_numpy())
        matched += len(positive)
        counts.append((pool, len(positive)))
    result.update(n_matched_reference_edges=matched, n_unmatched_reference_edges=excluded,
                  n_selected_prediction_edges=len(selected))
    if not matched:
        result["reason"] = "no_distance_observability_matched_unlabeled_background"
        return result, audit
    null = np.zeros(repeats, dtype=float)
    for pool, size in counts:
        recovery = np.array([edge in selected for edge in pool], dtype=int)
        for iteration in range(repeats):
            null[iteration] += recovery[rng.choice(len(pool), size=size, replace=False)].sum()
    expected = null.mean()
    result.update(status="evaluated", reason="", n_recovered_reference_edges=observed,
                  qualified_reference_recovery_fraction=observed / matched,
                  matched_unlabeled_recovery_fraction=expected / matched,
                  reference_recovery_enrichment=observed / expected if expected > 0 else "undefined_zero_background_recovery",
                  nominal_matched_background_tail_probability=(1 + int((null >= observed).sum())) / (repeats + 1),
                  interpretation="conditional_reference_recovery_against_matched_unlabeled_edges_not_false_positive_accuracy",
                  matched_background_definition=f"same_{anchor}_distance_bin_observability_bin_no_known_positive_edges")
    return result, audit


def write_table(path, values, default_columns):
    pd.DataFrame(values, columns=None if values else default_columns).to_csv(path, sep="\t", index=False)


def run(args):
    if args.output.exists():
        raise FileExistsError(f"Refusing overwrite: {args.output}")
    args.output.mkdir(parents=True, exist_ok=False)
    report = dict(status="not_run", reason="no_qualified_supplied_reference", scientific_adoption=False,
                  missing_edges="unlabeled_not_verified_negatives", seed=args.seed,
                  reference_types=sorted(TYPES), auroc_auprc="not_computed",
                  multiplicity="no_correction_nominal_descriptive_probabilities_no_scientific_calls",
                  supplied_reference_files=[str(path.resolve()) for path in args.reference])
    if not args.reference:
        (args.output / "status.json").write_text(json.dumps(report, indent=2) + "\n")
        write_table(args.output / "reference_recovery.tsv", [], GROUP + ["status", "reason"])
        return
    reference = validate_reference(pd.concat([table(path) for path in args.reference], ignore_index=True))
    if reference.duplicated(IDENTITY).any():
        raise ValueError("Reference duplicates: consolidate evidence per edge before evaluation")
    if reference.empty or args.universe is None:
        report["reason"] = "no_qualified_supplied_reference" if reference.empty else "qualified_candidate_universe_not_supplied"
        (args.output / "status.json").write_text(json.dumps(report, indent=2) + "\n")
        write_table(args.output / "reference_recovery.tsv", [], GROUP + ["status", "reason"])
        return
    universe = table(args.universe)
    require_columns(universe, IDENTITY + ["distance_bin", "observability_bin", "observable"], "Candidate universe")
    if universe.duplicated(IDENTITY).any():
        raise ValueError("Duplicate edge in candidate universe")
    if not set(universe.observable.str.lower()).issubset({"true", "false", "1", "0"}):
        raise ValueError("observable must be explicitly true/false or 1/0")
    if not set(universe.reference_type).issubset(TYPES):
        raise ValueError("Unknown universe reference type")
    priors = None if args.prior is None else table(args.prior)
    if priors is not None:
        require_columns(priors, IDENTITY + ["method", "prior_id"], "Prior manifest")
    link_paths = list(args.links)
    for network in args.network:
        if not (network / "completion.json").exists():
            raise ValueError("Network run has no completion marker")
        if json.loads((network / "completion.json").read_text()).get("status") != "completed" or any(network.glob("INVALIDATED*.json")):
            raise ValueError("Network run is incomplete or invalidated")
        link_paths.extend(sorted(network.glob("*/links.tsv")))
    predictions = [adapt_links(path) for path in link_paths if path.stat().st_size > 1]
    predictions = pd.concat(predictions, ignore_index=True) if predictions else pd.DataFrame(columns=GROUP + ["source_id", "target_id", "score"])
    if predictions.duplicated(GROUP + ["source_id", "target_id"]).any():
        raise ValueError("Duplicate predictions across supplied files")
    rng, results, audit = np.random.default_rng(args.seed), [], []
    for _, group in predictions.groupby(GROUP, sort=True):
        result, overlap = evaluate_group(group, reference, universe, priors, rng, args.top_k, args.background_repeats)
        results.append(result)
        audit.extend(overlap)
    observed_contexts = set(map(tuple, predictions[["reference_type", "cohort", "lineage"]].to_numpy()))
    for context, group in reference.groupby(["reference_type", "cohort", "lineage"], sort=True):
        if context not in observed_contexts:
            results.append(dict(zip(["reference_type", "cohort", "lineage"], context), method="", outer_fold="",
                status="not_run", reason="no_compatible_supplied_predictions", n_qualified_reference_edges=len(group)))
    write_table(args.output / "reference_recovery.tsv", results, GROUP + ["status", "reason"])
    write_table(args.output / "prior_overlap_audit.tsv", audit, GROUP + ["source_id", "target_id", "reference_id", "prior_overlap_state"])
    evaluated = sum(result["status"] == "evaluated" for result in results)
    report.update(status="completed" if results else "not_run", reason="" if results else "no_compatible_supplied_predictions",
                  n_evaluated_groups=evaluated, n_unevaluated_groups=len(results) - evaluated,
                  prior_audit="unassessed" if priors is None else "supplied_manifest",
                  tail_probability="nominal_descriptive_matched_background_probability_no_method_selection_or_causal_claims",
                  universe_path=str(args.universe.resolve()), prior_path="" if args.prior is None else str(args.prior.resolve()),
                  prediction_paths=[str(path.resolve()) for path in link_paths])
    (args.output / "status.json").write_text(json.dumps(report, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, action="append", default=[])
    parser.add_argument("--universe", type=Path)
    parser.add_argument("--prior", type=Path)
    parser.add_argument("--links", type=Path, action="append", default=[])
    parser.add_argument("--network", type=Path, action="append", default=[])
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--background-repeats", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.top_k < 1 or args.background_repeats < 1:
        parser.error("Positive --top-k and --background-repeats required")
    if args.reference and (not os.environ.get("SLURM_JOB_ID") or "login" in socket.gethostname()):
        raise RuntimeError("Reference evaluation requires a compute allocation")
    run(args)


if __name__ == "__main__":
    main()
