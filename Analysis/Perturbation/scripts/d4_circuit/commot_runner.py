#!/usr/bin/env python
"""COMMOT runner for D4 circuit. Subagent: circuit-runner-commot.

COMMOT propagates LR signaling on spatial data (Visium GSE192741). The runner
treats spatial neighborhoods as the receiver context.
"""
from __future__ import annotations
from _runner_template import D4Adapter, main_for_model


class CommotAdapter(D4Adapter):
    model_name = "commot"
    propagation_method = "commot"
    default_checkpoint = None  # COMMOT is not a checkpointed model

    def load_checkpoint(self) -> None:
        # TODO[circuit-runner-commot]: import commot; load Visium adata
        # (Analysis/Spatial/results/...) ; build LR graph.
        raise NotImplementedError("Fill load_checkpoint for COMMOT")

    def predict_for_hit(self, *, hit_row, context, substrate):
        # TODO[circuit-runner-commot]: knock out ligand, recompute signaling,
        # report receiver delta.
        raise NotImplementedError("Fill predict_for_hit for COMMOT")


if __name__ == "__main__":
    main_for_model(CommotAdapter)
