#!/usr/bin/env python3
"""Re-derive the sealed BG-001 graph inside the Resource candidate mirror."""

from __future__ import annotations

import json

from accepted_run_provenance import (
    EXPECTED_GRAPH_BYTES,
    EXPECTED_GRAPH_FILES,
    verify_accepted_graph,
)
from resource_contract import (
    ACCEPTED_RUN_INPUT,
    CANDIDATE_ID,
    ContractError,
    require,
)


def verify() -> dict[str, object]:
    graph = verify_accepted_graph(ACCEPTED_RUN_INPUT)
    require(len(graph) == EXPECTED_GRAPH_FILES, "mirrored graph count drift")
    total_bytes = sum(size for _, size, _ in graph.values())
    require(total_bytes == EXPECTED_GRAPH_BYTES, "mirrored graph byte drift")
    return {
        "status": "BG001_MIRRORED_PROVENANCE_VERIFIED",
        "candidate_id": CANDIDATE_ID,
        "graph_file_count": len(graph),
        "graph_size_bytes": total_bytes,
    }


if __name__ == "__main__":
    try:
        print(json.dumps(verify(), indent=2, sort_keys=True))
    except ContractError as error:
        raise SystemExit(f"MIRRORED PROVENANCE ERROR: {error}") from error
