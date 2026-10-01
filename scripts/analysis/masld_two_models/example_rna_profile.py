#!/usr/bin/env python3
"""Reproduce one supported raw RNA-to-profile call using a public source sample.

This example exercises the released interface. Its sample was in the model's
training source, so the output is not a held-participant accuracy estimate.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RELEASE = ROOT/"Analysis/MASLD_Model_Benchmark/release/masld-liver-chromatin-state-v1.2"
FIX = ROOT/"Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture"


def main(out):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Compute-node allocation required")
    out.mkdir(parents=True, exist_ok=False)
    participant = str(pd.read_csv(FIX/"molecular/participant_axis.tsv", sep="\t").participant_id.iloc[0])
    genes = pd.read_csv(FIX/"molecular/rna_feature_axis.tsv", sep="\t").stable_gene_id.astype(str).to_numpy()
    values = np.asarray(np.load(FIX/"molecular/rna_values.npy", mmap_mode="r")[0], float)
    table = pd.DataFrame({"stable_gene_id": genes, participant: values})
    counts = out/"example_fractional_RNA.tsv.gz"
    table.to_csv(counts, sep="\t", index=False, float_format="%.17g", compression="gzip")
    command = [sys.executable, str(RELEASE/"score.py"), "--counts", str(counts),
               "--out", str(out/"supporting_states.tsv"), "--profile", str(out/"predicted_H3K27ac.npz"),
               "--profile-form", "rrr_offset_cis", "--transport", "raw", "--json", str(out/"scorer_report.json")]
    with (out/"scorer.log").open("w") as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    module_spec = importlib.util.spec_from_file_location("released_rna_profile", RELEASE/"score.py")
    scorer = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(scorer)
    direct = scorer.profile(values[:, None], genes, form="rrr_offset_cis", transport="raw")
    with np.load(out/"predicted_H3K27ac.npz", allow_pickle=False) as cli:
        maximum = float(np.max(np.abs(cli["profile"]-direct["profile"])))
        if maximum > 1e-10 or cli["profile"].shape != (1, 96460) or float(cli["axis_coverage"]) != 1:
            raise ValueError("CLI example did not reproduce the direct released profile")
        if not np.array_equal(cli["region_key"].astype(str), direct["region_key"]):
            raise ValueError("CLI and direct region identities differ")
    report = {"sample_id": participant, "source": "GSE267145; already present in model training",
              "input_units": "deposited continuous RNA abundance estimates; no rounding",
              "input_orientation": "stable gene IDs by samples", "axis_coverage": 1.0,
              "output": "predicted concentration-residual H3K27ac log2 CPM on 96460 supported GRCh38 regions",
              "transport": "raw; fixed individual-sample transformation",
              "accuracy_claim": "none; this verifies the supported invocation and numerical reproduction",
              "CLI_vs_direct_max_absolute_difference": maximum,
              "command": command, "python": sys.version, "numpy": np.__version__,
              "slurm_job_id": os.environ["SLURM_JOB_ID"],
              "sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in (counts, RELEASE/"score.py", RELEASE/"weights/chromatin_state_v1_2.npz")}}
    (out/"example.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"sample": participant, "profile_shape": [1, 96460], "coverage": 1.0,
                      "CLI_vs_direct_max_absolute_difference": maximum}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args().out)
