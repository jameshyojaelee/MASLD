#!/usr/bin/env python3
"""Run a network-free synthetic GEARS core forward/backward fixture."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
import socket

import torch

from gears.model import GEARS_Model


class SyntheticBatch:
    def __init__(self, x: torch.Tensor, pert_idx: list[list[int]], batch: torch.Tensor):
        self.x = x
        self.pert_idx = pert_idx
        self.batch = batch


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def graph_edges(size: int) -> tuple[torch.Tensor, torch.Tensor]:
    sources = list(range(size)) + list(range(size - 1)) + list(range(1, size))
    targets = list(range(size)) + list(range(1, size)) + list(range(size - 1))
    edge_index = torch.tensor([sources, targets], dtype=torch.long)
    edge_weight = torch.ones(len(sources), dtype=torch.float32)
    return edge_index, edge_weight


def make_model(seed: int) -> GEARS_Model:
    torch.manual_seed(seed)
    gene_edges, gene_weights = graph_edges(8)
    pert_edges, pert_weights = graph_edges(3)
    return GEARS_Model(
        {
            "num_genes": 8,
            "num_perts": 3,
            "hidden_size": 8,
            "uncertainty": False,
            "num_go_gnn_layers": 1,
            "decoder_hidden_size": 4,
            "num_gene_gnn_layers": 1,
            "no_perturb": False,
            "device": "cpu",
            "G_coexpress": gene_edges,
            "G_coexpress_weight": gene_weights,
            "G_go": pert_edges,
            "G_go_weight": pert_weights,
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise RuntimeError("refusing to overwrite GEARS core fixture")
    torch.set_num_threads(1)
    torch.use_deterministic_algorithms(True)
    seed = 20260824
    torch.manual_seed(seed)
    x = torch.linspace(0.1, 1.6, 16, dtype=torch.float32).reshape(-1, 1)
    batch = SyntheticBatch(x, [[0], [1]], torch.tensor([0] * 8 + [1] * 8))
    original_socket = socket.socket
    socket.socket = lambda *unused_args, **unused_kwargs: (_ for _ in ()).throw(  # type: ignore[assignment]
        RuntimeError("runtime network is prohibited")
    )
    try:
        first = make_model(seed)
        first.eval()
        output = first(batch)
        if output.shape != (2, 8) or not torch.isfinite(output).all():
            raise RuntimeError("GEARS forward shape or finite-value contract failed")
        second = make_model(seed)
        second.eval()
        repeated = second(batch)
        if not torch.equal(output, repeated):
            raise RuntimeError("GEARS fixed-seed CPU forward is not deterministic")
        first.train()
        train_output = first(batch)
        loss = torch.mean((train_output - x.reshape(2, 8)) ** 2)
        loss.backward()
        gradients = [parameter.grad for parameter in first.parameters() if parameter.requires_grad]
        finite_gradients = [
            gradient for gradient in gradients if gradient is not None and torch.isfinite(gradient).all()
        ]
        if not finite_gradients or len(finite_gradients) != sum(value is not None for value in gradients):
            raise RuntimeError("GEARS backward finite-gradient contract failed")
    finally:
        socket.socket = original_socket
    source = Path(__import__("gears.model", fromlist=["unused"]).__file__).resolve()
    result = {
        "schema_version": "masld-bench-gears-core-forward-fixture-v1",
        "seed": seed,
        "device": "cpu",
        "torch_version": torch.__version__,
        "gears_model_source": source.as_posix(),
        "gears_model_source_sha256": digest(source),
        "output_shape": list(output.shape),
        "output_sum": float(output.sum().item()),
        "loss": float(loss.item()),
        "all_finite": bool(torch.isfinite(output).all()),
        "deterministic_exact_repeat": bool(torch.equal(output, repeated)),
        "finite_gradient_tensor_count": len(finite_gradients),
        "runtime_network_blocked": True,
        "external_GO_graph_used": False,
        "pickle_or_checkpoint_loaded": False,
        "biological_data_used": False,
        "scientific_admission": "blocked; core tensor fixture only",
        "gse313774_accessed": False,
    }
    if not math.isfinite(result["output_sum"]) or not math.isfinite(result["loss"]):
        raise RuntimeError("GEARS numeric receipt is not finite")
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "gears_core_forward.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
