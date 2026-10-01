#!/usr/bin/env python3
"""MANIFEST.tsv for the AD.1 gate run: every input, code file and output with its sha256.

Includes the local checkpoint files and the project FASTA, because the gate's whole claim is that a
particular set of weights on a particular reference reproduces a particular archived column.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import i1_common as C  # noqa: E402

CODE = [
    "i1_common.py",
    "i1_00_build_gate_set.py",
    "i1_00_build.sbatch",
    "i1_01_score_local.py",
    "i1_02_pilot.sbatch",
    "i1_03_pilot_analysis.py",
    "i1_04_gate_stats.py",
    "i1_05_zeroshot.py",
    "i1_06_gate.sbatch",
    "i1_07_manifest.py",
    "i1_08_write_results.py",
    "i1_09_final.sbatch",
    "i1_interim.py",
    "i1_smoke.sbatch",
    "ADAPTATION_ARM_PRESPEC.md",
    "ADAPTATION_ARM_AMENDMENT_01.md",
    "IMPLEMENTATION_SPEC.md",
]

INPUTS = [
    C.TIER_A4,
    C.C2_LABELS,
    C.C2_OOF,
    C.C2_METRICS,
    C.C2 / "tables/head_summary.tsv",
    C.C2 / "inputs/label_contract.json",
    pathlib.Path(C.FASTA_PATH),
    C.PROJ / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-weights-20260915T103442Z/acquisition.json",
    C.PROJ / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-weights-20260915T103442Z/checkpoint_inventory.tsv",
    C.PROJ / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-runtime-probe-20260915T104500Z/RESULTS.md",
    C.PROJ / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-runtime-probe-20260915T104500Z/ADDENDUM.md",
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True)
    args = ap.parse_args()
    run = pathlib.Path(args.run)
    src = pathlib.Path(__file__).resolve().parent

    rows = [("role", "path", "bytes", "sha256")]

    for name in CODE:
        p = src / name
        if p.exists():
            rows.append(("code", str(p), str(p.stat().st_size), C.sha256_file(p)))

    for p in INPUTS:
        p = pathlib.Path(p)
        if p.exists():
            rows.append(("input", str(p), str(p.stat().st_size), C.sha256_file(p)))

    for p in sorted(C.CHECKPOINT.rglob("*")):
        if p.is_file() and ".cache" not in p.parts:
            rows.append(("weights", str(p), str(p.stat().st_size), C.sha256_file(p)))

    for sub in ("inputs", "raw", "tables", "env"):
        for p in sorted((run / sub).rglob("*")):
            if p.is_file():
                rows.append(("output", str(p), str(p.stat().st_size), C.sha256_file(p)))

    for p in sorted((run / "logs").glob("*")):
        if p.is_file():
            rows.append(("log", str(p), str(p.stat().st_size), C.sha256_file(p)))

    (run / "MANIFEST.tsv").write_text("\n".join("\t".join(r) for r in rows) + "\n")
    C.log(f"wrote {run / 'MANIFEST.tsv'} with {len(rows) - 1} entries")

    C.write_json(
        run / "seeds.json",
        {
            "bootstrap_seed": C.BOOTSTRAP_SEED,
            "bootstrap_resamples": C.RESAMPLES,
            "comparison_seeds": list(C.COMPARISON_SEEDS),
            "comparison_seeds_use": (
                "the five comparison seeds index the fitted heads of the c2 package that this arm is "
                "paired against; the zero-shot arm has no fitted parameters, so no seed enters the "
                "local score and none of the five is consumed here"
            ),
            "resampling_unit": "1-Mb block (chr:floor(pos/1e6)), as deposited by the c2 label build",
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
