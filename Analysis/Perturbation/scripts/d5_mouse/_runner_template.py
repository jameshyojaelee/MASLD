"""D5 mouse runner template.

D5 task: re-run D1 / D2 mechanism prediction on the MOUSE atlas using mouse-
trained model weights. Emit GenePrediction rows keyed on the mouse gene symbol.
Cross-species concordance computed downstream by `concordance_scorer.py`.

USAGE:
    python state_mouse_runner.py --modality zero_shot --context HFD__C57BL6__M
"""
from __future__ import annotations

import sys
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import pandas as pd  # noqa: E402

from runner_base import (  # noqa: E402
    ModelAdapter,
    RESULTS_ROOT,
    build_arg_parser,
    load_hits,
    run_arm,
)

ARM = "D5"
DEFAULT_HITS = "d5_mouse_orthologs.csv"


class D5Adapter(ModelAdapter):
    arm = ARM
    model_name = "OVERRIDE_ME"
    default_checkpoint: str | None = None

    def load_checkpoint(self) -> None:
        raise NotImplementedError(
            f"{self.model_name}: implement load_checkpoint()"
        )

    def predict_for_hit(
        self,
        *,
        hit_row: pd.Series,
        context: str,
        substrate,
    ) -> list[dict]:
        """TODO[d5-runner]: emit GenePrediction dicts for mouse_gene perturbation.

        Schema matches D1 (GenePrediction). `target_gene` = mouse symbol.
        """
        raise NotImplementedError(
            f"{self.model_name}: implement predict_for_hit()"
        )


def main_for_model(adapter_cls: type[D5Adapter]) -> None:
    parser = build_arg_parser(adapter_cls.model_name, ARM, DEFAULT_HITS)
    args = parser.parse_args()

    hits = load_hits(args.hits, max_hits=args.max_hits)
    print(f"[{adapter_cls.model_name}] loaded {len(hits)} mouse ortholog hits")

    adapter = adapter_cls(
        modality=args.modality,
        checkpoint=args.checkpoint,
        dry_run=args.dry_run,
    )
    run_arm(
        adapter=adapter,
        hits_df=hits,
        context=args.context,
        out_dir=Path(args.out_dir),
        notes="D5 cross-species mouse mechanism",
    )
