"""LSEC (sinusoidal endothelial) receiver-response modeler. Subagent: circuit-receiver-lsec."""
from __future__ import annotations


class ReceiverModeler:
    receiver_cell_type = "LSEC"

    def __init__(self):
        # TODO[circuit-receiver-lsec]: load LSEC scVI model
        self.model = None

    def respond(self, hep_state):
        # TODO[circuit-receiver-lsec]
        raise NotImplementedError("Fill ReceiverModeler.respond for LSEC")
