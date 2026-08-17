#!/usr/bin/env python3
import argparse
import json

from masld_cl.config import load_config
from masld_cl.gpu_tolerance import measure_gpu_tolerance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--policy", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    result = measure_gpu_tolerance(load_config(args.config), args.policy, args.output)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
