"""D3 synergy runner template.

Each per-model runner (gears_runner.py, state_pair_runner.py,
scgpt_mask_runner.py, tahoe_drug_runner.py) subclasses D3Adapter.

D3 task: predict combinatorial KO synergy (pair / triple / quad) — emit
SynergyPrediction rows with `sigma_above_additive`. Consensus rule: 2σ above
additive in ≥2 models.

USAGE:
    python gears_runner.py --modality zero_shot --tier tier1_pairs --context all
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

ARM = "D3"


class D3Adapter(ModelAdapter):
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
        """TODO[d3-runner]: return SynergyPrediction dicts.

        Schema fields:
          gene1, gene2, gene3, gene4 (None if absent),
          additive_baseline, observed_double_or_higher,
          synergy_magnitude, synergy_class,
          sigma_above_additive, k562_bias_confidence.

        Example:
            g1, g2 = hit_row["gene1"], hit_row["gene2"]
            a = self.model.predict_ko(g1, substrate)
            b = self.model.predict_ko(g2, substrate)
            ab = self.model.predict_double_ko(g1, g2, substrate)
            additive = a + b
            sigma = (ab - additive) / self.model.noise_std
            return [dict(
                gene1=g1, gene2=g2,
                additive_baseline=float(additive.mean()),
                observed_double_or_higher=float(ab.mean()),
                synergy_magnitude=float(abs(ab - additive).mean()),
                synergy_class="synergistic" if sigma > 0 else "antagonistic",
                sigma_above_additive=float(sigma),
                k562_bias_confidence=0.5,
            )]
        """
        raise NotImplementedError(
            f"{self.model_name}: implement predict_for_hit()"
        )


def main_for_model(adapter_cls: type[D3Adapter]) -> None:
    parser = build_arg_parser(adapter_cls.model_name, ARM, "d3_tier1_pairs.csv")
    parser.add_argument(
        "--tier",
        default="tier1_pairs",
        choices=[
            "tier1_pairs",
            "tier2_pairs",
            "tier3_triples",
            "tier4_quadruples",
            "tier5_adaptive",
        ],
        help="Combinatorial tier to score",
    )
    args = parser.parse_args()

    # Resolve hit file from tier if --hits not explicitly overridden to a path that exists.
    from combinatorial_planner import HITS_DIR

    tier_path = HITS_DIR / f"d3_{args.tier}.csv"
    hits_path = args.hits if Path(args.hits).exists() and "tier" in args.hits else str(tier_path)
    hits = load_hits(hits_path, max_hits=args.max_hits)
    print(f"[{adapter_cls.model_name}] tier={args.tier} loaded {len(hits)} hits")

    adapter = adapter_cls(
        modality=args.modality,
        checkpoint=args.checkpoint,
        dry_run=args.dry_run,
    )

    run_arm(
        adapter=adapter,
        hits_df=hits,
        context=f"{args.tier}__{args.context}",
        out_dir=Path(args.out_dir),
        notes=f"D3 synergy / tier={args.tier}",
    )
