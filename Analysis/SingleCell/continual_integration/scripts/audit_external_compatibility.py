#!/usr/bin/env python3
"""Audit target donors, labels, and genes in the pinned HLiCA asset."""

from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config, write_json_exclusive
from masld_cl.external_compatibility import audit_external_compatibility
from masld_cl.reference_assessment import load_policy, sha256_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--label-map", required=True)
    parser.add_argument("--external-h5ad", required=True)
    parser.add_argument("--prepared-h5ad", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    policy = load_policy(args.policy)
    result = audit_external_compatibility(
        args.external_h5ad, args.prepared_h5ad, policy, args.label_map
    )
    result.update({
        "config_sha256": config["_config_sha256"],
        "policy_sha256": sha256_path(args.policy),
        "label_map_sha256": sha256_path(args.label_map),
        "external_h5ad_sha256": sha256_path(args.external_h5ad),
        "prepared_h5ad_sha256": sha256_path(args.prepared_h5ad),
    })
    write_json_exclusive(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
