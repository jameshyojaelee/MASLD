"""Shared utilities for the perturbation campaign.

Used by every arm's subagents:
- load_hepatocyte_atlas    : 657k hepatocyte AnnData + sex×stage×subtype stratification
- load_bulk_5cohort        : 847-sample bulk for bulk-compat models
- backtest_framework       : positive-control recovery curves; ejection gate
- output_schema            : pydantic schema for model predictions + consensus
- consensus_aggregator     : 2-of-N consensus per arm
- atlas_writer             : S8 columns → multi_evidence_atlas.csv
"""

from . import (  # noqa: F401
    atlas_writer,
    backtest_framework,
    consensus_aggregator,
    load_bulk_5cohort,
    load_hepatocyte_atlas,
    output_schema,
)
