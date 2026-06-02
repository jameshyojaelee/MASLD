#!/usr/bin/env python
"""scGen runner for D2 reversal. Subagent: reversal-runner-scgen."""
from __future__ import annotations
from _runner_template import D2Adapter, main_for_model


class ScGenAdapter(D2Adapter):
    model_name = "scgen"
    default_checkpoint = (
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
        "Analysis/Perturbation/results/finetuned_checkpoints/scgen_hep/"
    )

    def load_checkpoint(self) -> None:
        # TODO[reversal-runner-scgen]: scgen.SCGEN.load(self.checkpoint)
        raise NotImplementedError("Fill load_checkpoint for scGen")

    def predict_for_hit(self, *, hit_row, context, substrate):
        # TODO[reversal-runner-scgen]: compute reversal score for hit_row.gene
        # against self._ref_df (Diseased->Healthy signature).
        raise NotImplementedError("Fill predict_for_hit for scGen")


if __name__ == "__main__":
    main_for_model(ScGenAdapter)
