"""Macrophage receiver-response modeler (for STATE 2-stage).

Subagent: circuit-receiver-mac.
Loaded by `state_2stage_runner.py` via `importlib`.
"""
from __future__ import annotations


class ReceiverModeler:
    receiver_cell_type = "Macrophage"

    def __init__(self):
        # TODO[circuit-receiver-mac]: load macrophage scVI / state model so
        # `respond(modified_hep_state)` returns a delta-expression vector.
        self.model = None

    def respond(self, hep_state):
        # TODO[circuit-receiver-mac]: predict macrophage response to modified
        # hepatocyte secretome. Return object with .top_genes(k) and .magnitude.
        raise NotImplementedError("Fill ReceiverModeler.respond for Macrophage")
