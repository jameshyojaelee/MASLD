#!/usr/bin/env python
"""CPA (Compositional Perturbation Autoencoder) runner for D2 reversal.

Subagent: reversal-runner-cpa.
"""
from __future__ import annotations
from _runner_template import D2Adapter, main_for_model


class CPAAdapter(D2Adapter):
    model_name = "cpa"
    default_checkpoint = (
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
        "Analysis/Perturbation/results/finetuned_checkpoints/cpa_hep/"
    )

    def load_checkpoint(self) -> None:
        # TODO[reversal-runner-cpa]: from cpa import CPA
        raise NotImplementedError("Fill load_checkpoint for CPA")

    def predict_for_hit(self, *, hit_row, context, substrate):
        # TODO[reversal-runner-cpa]
        raise NotImplementedError("Fill predict_for_hit for CPA")


if __name__ == "__main__":
    main_for_model(CPAAdapter)
