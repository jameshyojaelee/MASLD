#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.metric_bundle import lock_metric_bundle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--contract-lock", required=True)
    parser.add_argument("--prepared", required=True)
    parser.add_argument("--prepared-lock", required=True)
    parser.add_argument("--selection-lock", required=True)
    parser.add_argument("--refit-lock", required=True)
    parser.add_argument("--confirmation-matrix-lock", required=True)
    parser.add_argument("--metric", action="append", required=True)
    parser.add_argument("--execution-record", action="append", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = lock_metric_bundle(
        load_config(args.config), args.contract_lock, args.prepared, args.prepared_lock,
        args.selection_lock, args.refit_lock, args.confirmation_matrix_lock,
        args.metric, args.execution_record, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
