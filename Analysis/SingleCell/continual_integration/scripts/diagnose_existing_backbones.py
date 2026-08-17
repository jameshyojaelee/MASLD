#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.backbone_diagnostic import diagnose_backbones
from masld_cl.config import load_config


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--reference-embedding", required=True)
    parser.add_argument("--harmony-embedding", required=True)
    parser.add_argument("--embedding", action="append", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = diagnose_backbones(
        load_config(args.config), args.reference_embedding, args.harmony_embedding,
        args.embedding, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
