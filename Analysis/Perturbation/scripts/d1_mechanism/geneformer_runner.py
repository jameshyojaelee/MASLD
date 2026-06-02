#!/usr/bin/env python
"""Geneformer (V1 + V2) runner for D1 mechanism.

Subagent: mechanism-runner-geneformer.
"""
from __future__ import annotations

from _runner_template import D1Adapter, main_for_model

DEFAULT_CHECKPOINT = (
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
    "Analysis/Perturbation/results/finetuned_checkpoints/geneformer_v2/"
)


class GeneformerAdapter(D1Adapter):
    model_name = "geneformer"
    default_checkpoint = DEFAULT_CHECKPOINT

    def load_checkpoint(self) -> None:
        # TODO[mechanism-runner-geneformer]:
        # from geneformer import InSilicoPerturber
        # self.model = InSilicoPerturber(...)
        raise NotImplementedError("Fill load_checkpoint for Geneformer")

    def predict_for_hit(self, *, hit_row, context, substrate):
        # TODO[mechanism-runner-geneformer]: run InSilicoPerturber, return
        # ranked downstream genes.
        raise NotImplementedError("Fill predict_for_hit for Geneformer")


if __name__ == "__main__":
    main_for_model(GeneformerAdapter)
