"""Campaign-v7 entry point for the frozen five-seed Cobolt evaluator."""

from __future__ import annotations

from . import rna_atac_cobolt_5seed_development as evaluator_v1


SOURCE_CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v5"
CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v7"


def main() -> int:
    """Reuse the audited evaluator while changing only its campaign binding."""

    if evaluator_v1.CAMPAIGN_ID != SOURCE_CAMPAIGN_ID:
        raise RuntimeError("base five-seed evaluator campaign binding changed")
    evaluator_v1.CAMPAIGN_ID = CAMPAIGN_ID
    try:
        return evaluator_v1.main()
    finally:
        evaluator_v1.CAMPAIGN_ID = SOURCE_CAMPAIGN_ID


if __name__ == "__main__":
    raise SystemExit(main())
