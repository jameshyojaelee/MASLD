"""Model-environment adapter protocol.

Adapters create assay-native predictions and artifacts.  Benchmark metrics are
computed only by :mod:`masld_bench.evaluators`.
"""

from .base import (
    Adapter,
    AdapterAction,
    AdapterError,
    AdapterReceipt,
    SubprocessAdapter,
    expected_action_dataset_ids,
    withheld_input_roles_for_action,
)

__all__ = [
    "Adapter",
    "AdapterAction",
    "AdapterError",
    "AdapterReceipt",
    "SubprocessAdapter",
    "expected_action_dataset_ids",
    "withheld_input_roles_for_action",
]
