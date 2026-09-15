"""Campaign-v12 entry point for the frozen five-seed Cobolt evaluator.

Binds the audited evaluator to the B6K campaign. Nothing about the evaluation
changes: the same frozen evaluator, the same registered strata roster, the same
metrics. Only the campaign identity moves, exactly as v6, v7 and v8 did.

DEVELOPMENT EVALUATOR, NOT A PRIMARY ENDPOINT. The v12 plan carries
metadata.primary_endpoint_scoring_allowed = None because the task has no variant
capability record, and selection.py raises SelectionError unless that field is
True. Best-model selection is therefore structurally unavailable for this campaign,
not merely withheld by policy. Output here is descriptive.
"""

from __future__ import annotations

from . import rna_atac_cobolt_5seed_development as evaluator_v1


SOURCE_CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v5"
CAMPAIGN_ID = "v1-rna-atac-cobolt-5seed-smoke-v12"


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
