#!/usr/bin/env python3
"""Numerical input sensitivity of the released RNA-to-profile scorer.

The 99 development participants are used only to exercise the released input
route. This is no estimate of held-person or external prediction accuracy.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RELEASE = ROOT/"Analysis/MASLD_Model_Benchmark/release/masld-liver-chromatin-state-v1.2"
FIX = ROOT/"Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture"
SEED = 20260922


def rms(a, b):
    return float(np.sqrt(np.mean((a-b)**2)))


def main(out):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    if out.exists():
        raise FileExistsError(out)
    out.mkdir(parents=True)
    spec = importlib.util.spec_from_file_location("released_chromatin_score", RELEASE/"score.py")
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    counts = np.asarray(np.load(FIX/"molecular/rna_values.npy", mmap_mode="r")[:20], float).T
    genes = pd.read_csv(FIX/"molecular/rna_feature_axis.tsv", sep="\t").stable_gene_id.astype(str).to_numpy()
    if counts.shape != (len(genes), 20):
        raise ValueError("RNA gene axis differs")
    full10 = scorer.profile(counts[:, :10], genes, form="rrr_offset_cis", transport="raw")
    single = scorer.profile(counts[:, :1], genes, form="rrr_offset_cis", transport="raw")
    exact_single = float(np.max(np.abs(single["profile"][0]-full10["profile"][0])))
    if exact_single > 1e-8:
        raise ValueError("Raw individual profile changed with co-scored people")
    marginal10 = scorer.profile(counts[:, :10], genes, form="rrr_offset_cis", transport="marginal")
    marginal20 = scorer.profile(counts, genes, form="rrr_offset_cis", transport="marginal")
    rows = [{"scenario": "raw individual versus ten-person batch", "coverage": 1.0,
             "rms_profile_difference": rms(single["profile"][0], full10["profile"][0]),
             "max_absolute_difference": exact_single},
            {"scenario": "marginal ten-person versus twenty-person batch", "coverage": 1.0,
             "rms_profile_difference": rms(marginal10["profile"], marginal20["profile"][:10]),
             "max_absolute_difference": float(np.max(np.abs(marginal10["profile"]-marginal20["profile"][:10])))}]
    for factor in (.5, 2.0, 4.0):
        perturbed = scorer.profile(counts[:, :1]*factor, genes, form="rrr_offset_cis", transport="raw")
        rows.append({"scenario": f"raw single-person library multiplied by {factor:g}", "coverage": 1.0,
                     "rms_profile_difference": rms(single["profile"], perturbed["profile"]),
                     "max_absolute_difference": float(np.max(np.abs(single["profile"]-perturbed["profile"])))})
    rng = np.random.default_rng(SEED)
    abundance = counts[:, :10].mean(1)
    for fraction, method in ((.9, "random"), (.75, "random"), (.6, "random"),
                             (.9, "lowest_abundance"), (.9, "highest_abundance")):
        size = int(round(fraction*len(genes)))
        if method == "random":
            keep = np.sort(rng.choice(len(genes), size=size, replace=False))
        elif method == "lowest_abundance":
            keep = np.sort(np.argsort(abundance)[:size])
        else:
            keep = np.sort(np.argsort(abundance)[-size:])
        output = scorer.profile(counts[keep, :10], genes[keep], form="rrr_offset_cis", transport="raw")
        rows.append({"scenario": f"raw {method} genes retained at {fraction:.0%}",
                     "coverage": float(output["axis_coverage"]),
                     "rms_profile_difference": rms(full10["profile"], output["profile"]),
                     "max_absolute_difference": float(np.max(np.abs(full10["profile"]-output["profile"])))})
    pd.DataFrame(rows).to_csv(out/"input_stability.tsv", sep="\t", index=False)
    (out/"scope.json").write_text(json.dumps({"status": "numerical scorer sensitivity only",
        "participants_exercised": 20, "training_source_reused": "GSE267145",
        "evaluation_claim": "none; participants are from model training source",
        "input": "continuous deposited RNA estimates, no integer rounding",
        "form": "rrr_offset_cis", "seed": SEED, "slurm_job_id": os.environ["SLURM_JOB_ID"]}, indent=2)+"\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True, type=Path)
    main(p.parse_args().out)
