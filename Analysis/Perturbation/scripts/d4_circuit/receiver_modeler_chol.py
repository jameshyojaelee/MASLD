"""Cholangiocyte receiver-response modeler. Subagent: circuit-receiver-chol."""
from __future__ import annotations


class ReceiverModeler:
    receiver_cell_type = "Cholangiocyte"

    def __init__(self):
        # TODO[circuit-receiver-chol]: load cholangiocyte scVI model
        self.model = None

    def respond(self, hep_state):
        # TODO[circuit-receiver-chol]
        raise NotImplementedError("Fill ReceiverModeler.respond for Cholangiocyte")
