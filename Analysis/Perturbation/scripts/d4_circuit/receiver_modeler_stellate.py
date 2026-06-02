"""Stellate (HSC) receiver-response modeler. Subagent: circuit-receiver-stellate."""
from __future__ import annotations


class ReceiverModeler:
    receiver_cell_type = "Stellate"

    def __init__(self):
        # TODO[circuit-receiver-stellate]: load HSC scVI model
        self.model = None

    def respond(self, hep_state):
        # TODO[circuit-receiver-stellate]
        raise NotImplementedError("Fill ReceiverModeler.respond for Stellate")
