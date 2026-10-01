#!/usr/bin/env python3
"""Write the bounded week-one recipes without launching any workload."""
import argparse
import csv
import itertools
import json
from pathlib import Path


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    recipes = []
    for length, pooling, head in itertools.product((2048, 16384), ("variant", "symmetric", "target"), ("ridge", "shared_mlp64")):
        recipes.append({"id": f"frozen_{length}_{pooling}_{head}", "mode": "frozen", "length": length,
            "pooling": pooling, "head": head, "hidden": 64, "loss": "mse", "seed": 1103,
            "alpha": 1000.0, "steps": 500,
            "status": "predeclared_development_comparison_not_selected_finalist"})
    for rank, depth, pooling in itertools.product((4, 16), (3, 5), ("variant", "symmetric")):
        recipes.append({"id": f"adapter_r{rank}_last{depth}_{pooling}", "mode": "adapter",
            "length": 2048, "pooling": pooling, "rank": rank, "last_blocks": depth,
            "hidden": 64, "loss": "mse", "seed": 1103, "backbone_lr": 3e-4,
            "steps": 100, "status": "single_split_seed_feasibility"})
    for rate in (3e-6, 1e-5):
        recipes.append({"id": f"partial_qv_last3_lr{rate:g}", "mode": "partial", "length": 2048,
            "pooling": "variant", "last_blocks": 3, "hidden": 64, "loss": "mse", "seed": 1103,
            "backbone_lr": rate, "steps": 100, "trainable_scope": "full_query_value_matrices_final_three_attention_blocks",
            "status": "single_split_seed_feasibility"})
    for recipe in recipes:
        path = out/(recipe["id"]+".json")
        if path.exists():
            raise FileExistsError(path)
        path.write_text(json.dumps(recipe, indent=2)+"\n")
    fields = sorted(set().union(*(r.keys() for r in recipes)))
    with (out/"recipes.tsv").open("w") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
        writer.writeheader(); writer.writerows(recipes)
    print(f"Wrote {len(recipes)} recipes: 12 frozen, 8 adapter, 2 partial; no jobs submitted")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    main(parser.parse_args().out)
