#!/usr/bin/env python3
"""Five-seed adaptation sweep at one validation fold, with fold 0 closed throughout.

Every adapter number deposited so far comes from one seed on one validation fold,
so arm ordering is confounded with seed noise and with which chromosomes happened
to validate. This runs the same twelve recipes at five seeds, and a sibling job
runs each of the four admissible validation folds.

Fold 0 is never trained on and never validated on, at any point. The trainer's
default pairing would rotate it into training as soon as the validation fold
moved, so the validation fold is passed explicitly and the held fold is pinned
at 0. Fold 0 stays available as the single untouched confirmation fold.

Selection is not performed here. This deposits fits; the nested selection reads
them afterwards and never reads fold 0.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

SEEDS = (20260917, 20260918, 20260919, 20260920, 20260921)
HELD_FOLD = 0


def configs_in(folder: Path) -> list:
    found = (sorted(folder.glob("adapter_*.json")) + sorted(folder.glob("partial_*.json"))
             + [folder / f"frozen_2048_{p}_shared_mlp64.json" for p in ("variant", "symmetric")])
    if len(found) != 12 or not all(p.exists() for p in found):
        raise ValueError("Expected 8 adapters, 2 partial updates and 2 matched frozen controls")
    return found


def main(args):
    if args.validation_fold == HELD_FOLD:
        raise ValueError("Fold 0 is closed; it cannot be the validation fold")
    args.out.mkdir(parents=True, exist_ok=False)
    trainer = Path(__file__).with_name("model_train.py")
    configs = configs_in(args.configs)
    started = time.monotonic()
    records = []

    def execute(config, seed):
        name = f"{config.stem}__seed{seed}"
        out = args.out / name
        command = [sys.executable, str(trainer), "--config", str(config), "--out", str(out),
                   "--steps", str(args.steps), "--microbatch", "4", "--validation-limit", "0",
                   "--held-fold", str(HELD_FOLD), "--validation-fold", str(args.validation_fold),
                   "--seed", str(seed), "--max-hours", str(args.per_arm_hours)]
        tick = time.monotonic()
        result = subprocess.run(command, check=False)
        record = {"name": name, "recipe": config.stem, "seed": seed,
                  "validation_fold": args.validation_fold, "held_fold": HELD_FOLD,
                  "exit_code": result.returncode, "wall_seconds": time.monotonic() - tick,
                  "steps_requested": args.steps, "config": str(config), "output": str(out)}
        completion = out / "feasibility.json"
        if completion.exists():
            done = json.loads(completion.read_text())["completed_steps"]
            record["steps_completed"] = done
            record["status"] = ("completed_fixed_exposure" if done == args.steps and not result.returncode
                                else "incomplete_feasibility_only")
        else:
            record["status"] = "failed_no_completed_feasibility_record"
        split = out / "split.json"
        if split.exists():
            s = json.loads(split.read_text())
            record["training_folds"] = s.get("training_folds")
            if HELD_FOLD in (s.get("training_folds") or []) or s.get("validation_fold") == HELD_FOLD:
                raise RuntimeError(f"{name} exposed the closed fold; stopping the sweep")
        records.append(record)
        (args.out / "runs.json").write_text(json.dumps(records, indent=2) + "\n")
        return record

    first = execute(configs[0], SEEDS[0])
    if first["status"] != "completed_fixed_exposure":
        raise RuntimeError("First arm did not complete its fixed exposure; sweep not admitted")
    total = len(configs) * len(SEEDS)
    projected = first["wall_seconds"] * total * args.safety
    remaining = args.max_hours * 3600 - (time.monotonic() - started)
    projection = {"validation_fold": args.validation_fold, "arms": total,
                  "first_arm_seconds": first["wall_seconds"], "safety_multiplier": args.safety,
                  "projected_seconds": projected, "remaining_seconds": remaining,
                  "admitted": projected < remaining + first["wall_seconds"] * args.safety}
    (args.out / "projection.json").write_text(json.dumps(projection, indent=2) + "\n")
    print(json.dumps(projection), flush=True)
    if not projection["admitted"]:
        (args.out / "runs.json").write_text(json.dumps(records, indent=2) + "\n")
        return

    for seed in SEEDS:
        for config in configs:
            if seed == SEEDS[0] and config == configs[0]:
                continue
            left = args.max_hours * 3600 - (time.monotonic() - started)
            if left < first["wall_seconds"] * args.safety + 300:
                records.append({"name": f"{config.stem}__seed{seed}", "recipe": config.stem,
                                "seed": seed, "status": "deferred_remaining_allocation"})
                continue
            execute(config, seed)

    done = [r for r in records if r.get("status") == "completed_fixed_exposure"]
    summary = {"validation_fold": args.validation_fold, "held_fold_closed": HELD_FOLD,
               "seeds": list(SEEDS), "recipes": [c.stem for c in configs],
               "planned_fits": total, "completed_fits": len(done),
               "deferred": sum(1 for r in records if r.get("status") == "deferred_remaining_allocation"),
               "wall_seconds": time.monotonic() - started,
               "selection_performed": False,
               "interpretation": "fits only; nested selection and fold-0 confirmation are separate steps"}
    (args.out / "sweep_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--configs", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--validation-fold", type=int, required=True)
    parser.add_argument("--steps", type=int, default=5000)
    parser.add_argument("--per-arm-hours", type=float, default=1.0)
    parser.add_argument("--max-hours", type=float, default=11.0)
    parser.add_argument("--safety", type=float, default=1.5)
    main(parser.parse_args())
