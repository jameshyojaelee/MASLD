#!/usr/bin/env python
"""scGPT mask-mode runner for D3 synergy. Subagent: synergy-runner-scgpt."""
from __future__ import annotations
from _runner_template import D3Adapter, main_for_model


class ScGPTMaskAdapter(D3Adapter):
    model_name = "scgpt_mask"
    default_checkpoint = (
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
        "Analysis/Perturbation/results/finetuned_checkpoints/scgpt_whole_human/"
    )

    def load_checkpoint(self) -> None:
        # TODO[synergy-runner-scgpt]
        raise NotImplementedError("Fill load_checkpoint for scGPT mask")

    def predict_for_hit(self, *, hit_row, context, substrate):
        # TODO[synergy-runner-scgpt]: gene masking gives KO-like delta;
        # multi-mask gives combinatorial.
        raise NotImplementedError("Fill predict_for_hit for scGPT mask")


if __name__ == "__main__":
    main_for_model(ScGPTMaskAdapter)
