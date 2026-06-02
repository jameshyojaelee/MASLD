#!/usr/bin/env python
"""STATE 2-stage runner for D4 circuit. Subagent: circuit-runner-state."""
from __future__ import annotations
from _runner_template import D4Adapter, main_for_model
import importlib


class State2StageAdapter(D4Adapter):
    model_name = "state_2stage"
    propagation_method = "state_2stage"
    default_checkpoint = (
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
        "Analysis/Perturbation/results/finetuned_checkpoints/state_se600m/"
    )

    def load_checkpoint(self) -> None:
        # TODO[circuit-runner-state]: load STATE; also instantiate per-receiver
        # response modelers from receiver_modeler_{mac,stellate,lsec,chol}.
        self.receiver_modelers = {}
        for ct in ("mac", "stellate", "lsec", "chol"):
            try:
                mod = importlib.import_module(f"receiver_modeler_{ct}")
                self.receiver_modelers[ct] = mod.ReceiverModeler()
            except (ImportError, AttributeError):
                self.receiver_modelers[ct] = None
        raise NotImplementedError("Fill load_checkpoint for STATE 2-stage")

    def predict_for_hit(self, *, hit_row, context, substrate):
        # TODO[circuit-runner-state]
        raise NotImplementedError("Fill predict_for_hit for STATE 2-stage")


if __name__ == "__main__":
    main_for_model(State2StageAdapter)
