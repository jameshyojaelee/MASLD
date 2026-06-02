"""LOO-dataset Hotspot re-runs + module stability scoring.

For one (cell_type, dataset_holdout) pair: re-run 501's pipeline with that
dataset excluded; later, compare module assignments to the full-run via
503b_loo_consolidate.py.

Usage:
    python 503_loo_stability.py --cell-type hepatocytes --dataset-holdout GSE136103
"""
from __future__ import annotations
import argparse
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from hotspot_io import RESULTS, RUN_ORDER  # noqa: F401  (RESULTS used by consolidator)

# Reuse 501's machinery without re-implementing the pipeline
_spec = importlib.util.spec_from_file_location(
    "_h501", str(Path(__file__).parent / "501_run_hotspot.py")
)
_h501 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_h501)


def main(cell_type: str, dataset_holdout: str) -> None:
    print(f"[LOO] {cell_type}  hold-out: {dataset_holdout}")
    _h501.main(cell_type=cell_type, smoke=False, loo_dataset=dataset_holdout)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="LOO-dataset stability: re-run Hotspot with one dataset excluded."
    )
    ap.add_argument(
        "--cell-type",
        required=True,
        choices=RUN_ORDER,
        help="Cell-type subset to analyse.",
    )
    ap.add_argument(
        "--dataset-holdout",
        required=True,
        help=(
            "Dataset name to exclude (must match 'dataset' column in donor_metadata.tsv). "
            "Valid values: GSE136103 GSE174748 GSE185477 GSE189600 GSE202379 GSE244832 Liver_Atlas"
        ),
    )
    main(**vars(ap.parse_args()))
