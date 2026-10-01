#!/usr/bin/env python3
"""Serial, budget-bounded microbatch-4 feasibility and 8+2 adaptation comparison.

A fresh 100-step throughput probe includes loading, checks, host batching,
checkpointing and validation. A conservative projection must fit before the
eight adapter, two partial-q/v, and two matched frozen-head runs are admitted.
No protected outcomes or final model selection enter this development sweep.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    trainer = Path(__file__).with_name("model_train.py")
    records = []
    def execute(config, name, steps, validation_limit, max_hours):
        out = args.out/name
        command = [sys.executable, str(trainer), "--config", str(config), "--out", str(out),
            "--steps", str(steps), "--microbatch", "4", "--validation-limit", str(validation_limit),
            "--max-hours", str(max_hours)]
        tick = time.monotonic()
        result = subprocess.run(command, check=False)
        record = {"name": name, "exit_code": result.returncode, "wall_seconds": time.monotonic()-tick,
                  "steps_requested": steps, "config": str(config), "output": str(out)}
        completion = out/"feasibility.json"
        if completion.exists():
            actual = json.loads(completion.read_text())["completed_steps"]
            record["steps_completed"] = actual
            record["status"] = "completed_fixed_exposure" if actual == steps and not result.returncode else "incomplete_feasibility_only"
        else:
            record["status"] = "failed_no_completed_feasibility_record"
        records.append(record)
        (args.out/"runs.json").write_text(json.dumps(records, indent=2)+"\n")
        return record
    probe = execute(args.configs/"adapter_r4_last3_variant.json", "probe_microbatch4", 100, 32, min(1.0,args.max_hours))
    if probe["exit_code"] or probe.get("steps_completed") != 100:
        raise RuntimeError("Microbatch4 probe failed; full sweep not admitted")
    feasibility = json.loads((args.out/"probe_microbatch4/feasibility.json").read_text())
    import pandas as pd
    root = Path(__file__).resolve().parents[3]
    labels = pd.read_csv(root/"GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz", sep="\t", usecols=["heldout_fold"])
    valid_n = int((labels.heldout_fold == 1).sum())
    train_cost = feasibility["iteration_wall_seconds_median"]*args.steps
    evaluation_cost = (feasibility["validation_first_batch_seconds"] +
        feasibility["validation_warm_seconds_per_example"]*max(0,valid_n-4))
    overhead = max(0, probe["wall_seconds"]-feasibility["training_wall_seconds_including_batching_and_progress_writes"]-feasibility["validation_wall_seconds"])
    # Twice measured costs cover deeper/rank16 and full-q/v pilots and filesystem variation.
    per_arm = 2*(train_cost+evaluation_cost+overhead)
    required = 12*per_arm
    remaining = args.max_hours*3600-(time.monotonic()-started)
    projection = {"probe_total_seconds": probe["wall_seconds"], "validation_rows": valid_n,
        "per_arm_training_seconds": train_cost, "per_arm_evaluation_seconds": evaluation_cost,
        "per_arm_fixed_overhead_seconds": overhead, "safety_multiplier": 2,
        "projected_12_arm_seconds": required, "remaining_seconds": remaining,
        "admitted": required < remaining-300}
    (args.out/"projection.json").write_text(json.dumps(projection, indent=2)+"\n")
    print(json.dumps(projection), flush=True)
    if not projection["admitted"]:
        return
    configs = sorted(args.configs.glob("adapter_*.json"))+sorted(args.configs.glob("partial_*.json"))
    configs += [args.configs/f"frozen_2048_{p}_shared_mlp64.json" for p in ("variant", "symmetric")]
    if len(configs) != 12:
        raise ValueError("Expected 8 adapters, 2 partials and 2 matched frozen controls")
    for config in configs:
        remaining = args.max_hours*3600-(time.monotonic()-started)
        if remaining < per_arm+300:
            records.append({"name": config.stem, "status": "deferred_remaining_budget"})
            continue
        execute(config, config.stem, args.steps, 0, min(remaining-120,per_arm+600)/3600)
    (args.out/"runs.json").write_text(json.dumps(records, indent=2)+"\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--max-hours", type=float, default=5.8)
    main(parser.parse_args())
