#!/usr/bin/env python3
"""Describe effect-magnitude capture within already selected Currin leads.

This asks which already measured association leads a fixed model would choose
at a fixed budget. It does not measure functional-variant enrichment, assay
success, new-variant transfer, MASLD risk, or a target-gene link.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
DESIGN = Path(__file__).with_name("allele_prioritization_addendum.json")
LABELS = ROOT/"GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz"
MODELS = ("native1m", "adapter", "combined")
BUDGETS = (.1, .2, .3)
SEED = 20260923


def source_peak_metadata(path=LABELS):
    membership = pd.read_csv(path, sep="\t", usecols=["heldout_fold"])
    if not membership.heldout_fold.isin(range(5)).all():
        raise ValueError("Unexpected source fold")
    skip = np.flatnonzero(~membership.heldout_fold.isin(range(1, 5)).to_numpy())+1
    fields = ["lead_variant_id", "chr", "pos_hg38", "ref", "alt", "heldout_fold",
              "peak_id", "peak_start_hg38", "peak_stop_hg38"]
    source = pd.read_csv(path, sep="\t", skiprows=skip, usecols=fields)
    source["key"] = source.lead_variant_id.str.removeprefix("chr")
    if source.key.duplicated().any() or not source.heldout_fold.isin(range(1, 5)).all():
        raise ValueError("Source identity or fold differs")
    return source


def read_predictions(path):
    membership = pd.read_csv(path, sep="\t", usecols=["fold"])
    if not membership.fold.isin(range(1, 5)).all():
        raise ValueError("Prediction file contains a forbidden fold; effects not parsed")
    frame = pd.read_csv(path, sep="\t")
    if frame.key.duplicated().any() or frame.groupby("chr").fold.nunique().max() != 1:
        raise ValueError("Prediction identity or chromosome-fold mapping differs")
    if not np.isfinite(frame[["beta_alt", *MODELS]].to_numpy(float)).all():
        raise ValueError("Nonfinite measured effect or prediction")
    return frame


def ranked_selection(frame):
    """Selection uses only identity, fold, and predictions, never measured beta."""
    result = frame.copy()
    result["selection_tie_hash"] = result.key.map(lambda key: hashlib.sha256(
        f"{SEED}|{key}".encode()).hexdigest())
    for model in MODELS:
        rank = pd.Series(index=result.index, dtype=np.int64)
        for _, part in result.groupby("fold", sort=True):
            order = part.assign(absolute_prediction=part[model].abs()).sort_values(
                ["absolute_prediction", "selection_tie_hash"], ascending=[False, True]).index
            rank.loc[order] = np.arange(1, len(order)+1)
        result[model+"_absolute_rank_within_fold"] = rank.astype(int)
        for budget in BUDGETS:
            limit = result.groupby("fold").key.transform("size").map(lambda n: int(np.ceil(n*budget)))
            result[f"{model}_selected_{int(budget*100)}pct"] = rank <= limit
    return result


def contributions(frame):
    records = []
    for model in MODELS:
        for budget in BUDGETS:
            for chrom, part in frame.groupby("chr", sort=True):
                selected = part[f"{model}_selected_{int(budget*100)}pct"].to_numpy(bool)
                y = part.beta_alt.to_numpy(float)
                pred = part[model].to_numpy(float)
                nonzero = y != 0
                correct = (np.sign(y) == np.sign(pred)) & nonzero
                fold_n = int(frame.fold.eq(part.fold.iloc[0]).sum())
                probability = np.ceil(budget*fold_n)/fold_n
                records.append({"model": model, "budget_fraction": budget, "chr": chrom,
                                "fold": int(part.fold.iloc[0]), "population_n": len(part),
                                "total_abs_beta": np.abs(y).sum(), "selected_n": int(selected.sum()),
                                "selected_abs_beta": np.abs(y[selected]).sum(),
                                "selected_direction_correct": int((selected & correct).sum()),
                                "selected_direction_n": int((selected & nonzero).sum()),
                                "random_expected_selected_n": probability*len(part),
                                "random_expected_abs_beta": probability*np.abs(y).sum(),
                                "random_expected_direction_correct": probability*correct.sum(),
                                "random_expected_direction_n": probability*nonzero.sum()})
    return pd.DataFrame(records)


def estimates(part, weights):
    def total(column):
        return weights@part[column].to_numpy(float)
    capture = total("selected_abs_beta")/total("total_abs_beta")
    random_capture = total("random_expected_abs_beta")/total("total_abs_beta")
    direction = np.divide(total("selected_direction_correct"), total("selected_direction_n"),
                          out=np.full(weights.shape[0], np.nan), where=total("selected_direction_n") > 0)
    random_direction = np.divide(total("random_expected_direction_correct"), total("random_expected_direction_n"),
                                 out=np.full(weights.shape[0], np.nan), where=total("random_expected_direction_n") > 0)
    return {"captured_abs_beta_fraction": capture,
            "uniform_random_expected_capture": random_capture,
            "capture_minus_uniform_random": capture-random_capture,
            "selected_direction_concordance": direction,
            "same_model_uniform_random_expected_direction": random_direction,
            "direction_minus_same_model_uniform_random": direction-random_direction}


def summarize(grouped, draws=4000):
    chromosomes = sorted(grouped.chr.unique())
    rng = np.random.default_rng(SEED)
    # Identical cluster multiplicities across all models/budgets make every
    # reported difference paired. Selection flags are held fixed in all draws.
    weights = rng.multinomial(len(chromosomes), np.repeat(1/len(chromosomes), len(chromosomes)), size=draws)
    rows = []
    for (model, budget), part in grouped.groupby(["model", "budget_fraction"], sort=True):
        part = part.set_index("chr").loc[chromosomes].reset_index()
        point = estimates(part, np.ones((1, len(chromosomes))))
        sampled = estimates(part, weights)
        row = {"model": model, "budget_fraction": budget, "selected_n": int(part.selected_n.sum()),
               "population_n": int(part.population_n.sum()), "n_chromosomes": len(chromosomes),
               "actual_selected_fraction": part.selected_n.sum()/part.population_n.sum(),
               "selected_direction_n": int(part.selected_direction_n.sum())}
        for metric in point:
            row[metric] = float(point[metric][0])
            finite = sampled[metric][np.isfinite(sampled[metric])]
            lower, upper = np.quantile(finite, [.025, .975]) if len(finite) else (np.nan, np.nan)
            row[metric+"_pointwise_CI95_lower"] = float(lower)
            row[metric+"_pointwise_CI95_upper"] = float(upper)
        rows.append(row)
    return pd.DataFrame(rows)


def main(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    if args.out.exists():
        raise FileExistsError(args.out)
    design = json.loads(DESIGN.read_text())
    if design["models"] != list(MODELS) or design["budgets"] != list(BUDGETS) or design["bootstrap_seed"] != SEED:
        raise ValueError("Ranking implementation differs from fixed addendum")
    frame = read_predictions(args.predictions)
    if len(frame) != 24430 or set(frame.fold) != set(range(1, 5)) or frame.chr.nunique() != 17:
        raise ValueError("Matched population differs from the fixed development census")
    metadata = source_peak_metadata()
    if len(metadata) != 25748:
        raise ValueError("Original development population differs")
    joined = frame.merge(metadata, on="key", suffixes=("", "_source"), validate="one_to_one")
    if len(joined) != len(frame):
        raise ValueError("Source peak metadata coverage is incomplete")
    for field in ("chr", "ref", "alt"):
        if not joined[field].equals(joined[field+"_source"]):
            raise ValueError("Source allele identity differs: "+field)
    if not np.array_equal(joined.fold, joined.heldout_fold):
        raise ValueError("Source fold identity differs")
    ranked = ranked_selection(joined)
    grouped = contributions(ranked)
    summary = summarize(grouped)
    args.out.mkdir(parents=True)
    fields = ["key", "chr", "pos_hg38", "ref", "alt", "fold", "peak_id",
              "peak_start_hg38", "peak_stop_hg38", "beta_alt", *MODELS]
    fields += [column for column in ranked if "_absolute_rank_within_fold" in column or "_selected_" in column]
    table = ranked[fields].copy()
    table["target_gene_evidence"] = "unsupported"
    table["prediction_state"] = "inspected Currin development; significance-selected leads"
    table["proposed_experiment"] = "Test the allele substitution at the variant and measure accessibility of the source peak; assess target-gene links separately"
    table.to_csv(args.out/"prioritization.tsv.gz", sep="\t", index=False)
    selected20 = table[[f"{model}_selected_20pct" for model in MODELS]].any(axis=1)
    table.loc[selected20].to_csv(args.out/"primary_20pct_selection_union.tsv.gz", sep="\t", index=False)
    grouped.to_csv(args.out/"fixed_selection_chromosome_contributions.tsv", sep="\t", index=False)
    summary.to_csv(args.out/"selection_summary.tsv", sep="\t", index=False)
    result = {"evidence_state": "secondary descriptive development use test; no superiority inference",
              "population": "significance-selected Currin association leads",
              "matched_n": len(table), "original_development_n": len(metadata),
              "coverage": len(table)/len(metadata), "chromosomes": frame.chr.nunique(),
              "primary_budget": .2, "primary_20pct": summary.loc[summary.budget_fraction.eq(.2)].to_dict("records"),
              "uncertainty": design["uncertainty"], "intervals": "pointwise 95%; not simultaneous",
              "hypothesis_tests": "none", "fold0_effects_parsed": 0,
              "target_gene_evidence": "unsupported", "limits": design["known_limits"],
              "hashes": {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                         for path in (DESIGN, Path(__file__), args.predictions)},
              "environment": {"python": sys.version, "numpy": np.__version__, "pandas": pd.__version__,
                              "host": platform.node(), "slurm_job_id": os.environ["SLURM_JOB_ID"]}}
    (args.out/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    (args.out/"design.json").write_text(json.dumps(design, indent=2)+"\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args())
