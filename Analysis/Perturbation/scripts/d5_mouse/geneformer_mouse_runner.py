#!/usr/bin/env python
"""Geneformer (mouse) runner for D5. Subagent: mouse-runner-geneformer."""
from __future__ import annotations
from _runner_template import D5Adapter, main_for_model


class GeneformerMouseAdapter(D5Adapter):
    model_name = "geneformer_mouse"
    default_checkpoint = (
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
        "Analysis/Perturbation/results/finetuned_checkpoints/geneformer_mouse/"
    )

    def load_checkpoint(self) -> None:
        # TODO[mouse-runner-geneformer]: Geneformer ships a mouse model variant
        raise NotImplementedError("Fill load_checkpoint for Geneformer mouse")

    def predict_for_hit(self, *, hit_row, context, substrate):
        # TODO[mouse-runner-geneformer]
        raise NotImplementedError("Fill predict_for_hit for Geneformer mouse")


if __name__ == "__main__":
    main_for_model(GeneformerMouseAdapter)
