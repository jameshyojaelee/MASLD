#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.training import train_update


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--prepared-lock", required=True)
    parser.add_argument("--contract-lock", required=True)
    parser.add_argument("--reference-model", required=True)
    parser.add_argument("--reference-embedding")
    parser.add_argument("--policy")
    parser.add_argument("--output", required=True)
    parser.add_argument("--model-kind", required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--ewc-lambda", type=float, required=True)
    parser.add_argument("--distillation-weight", type=float, default=0.0)
    parser.add_argument("--replay-fraction", type=float, required=True)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--production", action="store_true")
    parser.add_argument("--selection-lock")
    parser.add_argument("--refit-lock")
    parser.add_argument(
        "--replay-mode", default="random",
        choices=("random", "bi_bottom", "bi_top", "bi_step"),
    )
    parser.add_argument("--bi-replay-unlock")
    parser.add_argument("--query-dataset", action="append")
    parser.add_argument("--control-fisher-dataset", action="append")
    parser.add_argument("--held-out-dataset", action="append")
    parser.add_argument("--no-export-embedding", action="store_true")
    args = parser.parse_args()
    result = train_update(
        load_config(args.config), args.prepared, args.contract_lock, args.prepared_lock, args.reference_model,
        args.output, args.model_kind, args.method, args.ewc_lambda,
        args.replay_fraction, args.seed, args.production,
        selection_lock=args.selection_lock,
        query_datasets=args.query_dataset,
        control_fisher_datasets=args.control_fisher_dataset,
        held_out_datasets=args.held_out_dataset,
        export_embedding=not args.no_export_embedding,
        refit_lock=args.refit_lock,
        replay_mode=args.replay_mode,
        bi_replay_unlock=args.bi_replay_unlock,
        distillation_weight=args.distillation_weight,
        reference_embedding_path=args.reference_embedding,
        policy_value=args.policy,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
