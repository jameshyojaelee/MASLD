#!/usr/bin/env python
"""scGen mouse runner for D5. Subagent: mouse-runner-scgen."""
from __future__ import annotations
from _runner_template import D5Adapter, main_for_model


class ScGenMouseAdapter(D5Adapter):
    model_name = "scgen_mouse"
    default_checkpoint = (
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
        "Analysis/Perturbation/results/finetuned_checkpoints/scgen_mouse/"
    )

    def load_checkpoint(self) -> None:
        # TODO[mouse-runner-scgen]
        raise NotImplementedError("Fill load_checkpoint for scGen mouse")

    def predict_for_hit(self, *, hit_row, context, substrate):
        # TODO[mouse-runner-scgen]
        raise NotImplementedError("Fill predict_for_hit for scGen mouse")


if __name__ == "__main__":
    main_for_model(ScGenMouseAdapter)
