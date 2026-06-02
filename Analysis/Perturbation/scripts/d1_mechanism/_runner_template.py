"""D1 mechanism runner — per-model template.

Each per-model runner (state_runner.py, scgpt_runner.py, geneformer_runner.py,
tahoe_runner.py) imports this template and customizes:

  - MODEL_NAME
  - DEFAULT_CHECKPOINT
  - load_checkpoint(self)
  - predict_for_hit(self, hit_row, context, substrate)

D1 task: for each hit, return top-100 downstream genes with predicted logFC.

USAGE:
    python state_runner.py --modality zero_shot --context F__Steatohepatitis__Progressor
    python state_runner.py --modality fine_tuned --context F__Steatohepatitis__Progressor

Output: results/d1_mechanism/<model>_<modality>_<context>.json (ModelRunOutput JSON).
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

ARM = "D1"
DEFAULT_HITS = "d1_mechanism_hits.csv"
TOP_N_DOWNSTREAM = 100


class D1Adapter(ModelAdapter):
    """Mechanism prediction — emit GenePrediction rows per hit."""

    arm = ARM
    model_name = "OVERRIDE_ME"
    default_checkpoint: str | None = None

    # ----- TODO: subagent fills these in ---------------------------------
    def load_checkpoint(self) -> None:
        """TODO[mechanism-runner]: load weights from self.checkpoint.

        Example (STATE):
            from state_sdk import StateModel
            self.model = StateModel.from_pretrained(self.checkpoint)
        """
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
        """TODO[mechanism-runner]: call self.model on (hit_row.gene, substrate)
        and return TOP_N_DOWNSTREAM GenePrediction dicts.

        Each dict must match `output_schema.GenePrediction`:
          target_gene, downstream_gene, logFC_predicted, abs_rank, direction, confidence.

        Example:
            target = hit_row["gene"]
            ranks = self.model.predict_downstream(target, substrate=substrate, k=TOP_N_DOWNSTREAM)
            return [
                dict(
                    target_gene=target,
                    downstream_gene=g.symbol,
                    logFC_predicted=g.logfc,
                    abs_rank=i + 1,
                    direction="up" if g.logfc > 0 else "down",
                    confidence=g.score,
                )
                for i, g in enumerate(ranks)
            ]
        """
        raise NotImplementedError(
            f"{self.model_name}: implement predict_for_hit()"
        )


def main_for_model(adapter_cls: type[D1Adapter]) -> None:
    parser = build_arg_parser(adapter_cls.model_name, ARM, DEFAULT_HITS)
    args = parser.parse_args()

    hits = load_hits(args.hits, max_hits=args.max_hits)
    print(f"[{adapter_cls.model_name}] loaded {len(hits)} hits from {args.hits}")

    adapter = adapter_cls(
        modality=args.modality,
        checkpoint=args.checkpoint,
        dry_run=args.dry_run,
    )

    out_dir = Path(args.out_dir)
    run_arm(
        adapter=adapter,
        hits_df=hits,
        context=args.context,
        out_dir=out_dir,
        notes=f"D1 mechanism / top-{TOP_N_DOWNSTREAM} downstream",
    )
