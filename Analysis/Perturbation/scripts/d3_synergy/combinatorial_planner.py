"""D3 combinatorial search planner.

Decides which tier of combinatorial hits to score given a compute budget.

Tiers (from data/hits/):
  - d3_tier1_pairs.csv      : high-priority pairs (~10k)
  - d3_tier2_pairs.csv      : extended pairs    (~240k)
  - d3_tier3_triples.csv    : seed-driven triples
  - d3_tier4_quadruples.csv : seed-driven quads
  - d3_tier5_adaptive.csv   : adaptive expansion (filled after tier 1 results)

Usage:
    from combinatorial_planner import plan_run
    plan = plan_run(budget_hours=24, gpu_hours_per_pair=0.001)
    # -> list of (tier_name, hit_path, max_rows) to feed to runners.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
HITS_DIR = PROJECT_ROOT / "Analysis/Perturbation/data/hits"

TIERS = [
    ("tier1_pairs", HITS_DIR / "d3_tier1_pairs.csv"),
    ("tier2_pairs", HITS_DIR / "d3_tier2_pairs.csv"),
    ("tier3_triples", HITS_DIR / "d3_tier3_triples.csv"),
    ("tier4_quadruples", HITS_DIR / "d3_tier4_quadruples.csv"),
    ("tier5_adaptive", HITS_DIR / "d3_tier5_adaptive.csv"),
]


@dataclass
class PlanItem:
    tier_name: str
    hit_path: Path
    max_rows: int
    estimated_gpu_hours: float


def plan_run(
    *,
    budget_hours: float = 24,
    gpu_hours_per_pair: float = 1.0 / 1000,
    triple_multiplier: float = 3.0,
    quadruple_multiplier: float = 6.0,
) -> list[PlanItem]:
    """Greedy: fill tier1 first, then tier2, then triples, then quads.

    Returns one PlanItem per tier with max_rows capped by remaining budget.
    """
    import pandas as pd

    remaining = budget_hours
    plan: list[PlanItem] = []
    for name, path in TIERS:
        if not path.exists():
            continue
        n_rows = len(pd.read_csv(path))
        if "triples" in name:
            cost = gpu_hours_per_pair * triple_multiplier
        elif "quad" in name:
            cost = gpu_hours_per_pair * quadruple_multiplier
        else:
            cost = gpu_hours_per_pair
        max_rows = int(min(n_rows, remaining / max(cost, 1e-6)))
        if max_rows <= 0:
            continue
        plan.append(
            PlanItem(
                tier_name=name,
                hit_path=path,
                max_rows=max_rows,
                estimated_gpu_hours=max_rows * cost,
            )
        )
        remaining -= max_rows * cost
        if remaining <= 0:
            break
    return plan


if __name__ == "__main__":
    for p in plan_run():
        print(p)
