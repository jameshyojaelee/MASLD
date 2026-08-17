#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from masld_cl.config import load_config
from masld_cl.final_adoption_decision import write_final_adoption_decision


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--test-log", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = write_final_adoption_decision(
        load_config(args.config), args.policy, args.test_log, args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
