#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.external_reference_common import prepare_external26_common


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--contract-lock", required=True)
    parser.add_argument("--canonical-prepared", required=True)
    parser.add_argument("--canonical-prepared-lock", required=True)
    parser.add_argument("--external-h5ad", required=True)
    parser.add_argument("--external-compatibility", required=True)
    parser.add_argument("--label-map", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = prepare_external26_common(
        load_config(args.config), args.contract_lock, args.canonical_prepared,
        args.canonical_prepared_lock, args.external_h5ad, args.external_compatibility,
        args.label_map, args.policy, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
