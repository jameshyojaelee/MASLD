#!/usr/bin/env python
"""scGPT (whole-human checkpoint) runner for D1 mechanism.

Subagent: mechanism-runner-scgpt.
"""
from __future__ import annotations

from _runner_template import D1Adapter, main_for_model

DEFAULT_CHECKPOINT = (
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
    "Analysis/Perturbation/results/finetuned_checkpoints/scgpt_whole_human/"
)


class ScGPTAdapter(D1Adapter):
    model_name = "scgpt"
    default_checkpoint = DEFAULT_CHECKPOINT

    def load_checkpoint(self) -> None:
        # TODO[mechanism-runner-scgpt]: load scGPT (e.g. `scgpt.tasks.GeneEmbedder`)
        raise NotImplementedError("Fill load_checkpoint for scGPT")

    def predict_for_hit(self, *, hit_row, context, substrate):
        # TODO[mechanism-runner-scgpt]: use scGPT perturbation prediction head
        raise NotImplementedError("Fill predict_for_hit for scGPT")


if __name__ == "__main__":
    main_for_model(ScGPTAdapter)
