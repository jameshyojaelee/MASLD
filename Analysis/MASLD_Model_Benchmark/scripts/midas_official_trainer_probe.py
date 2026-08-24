#!/usr/bin/env python3
"""Probe the exact scmidas 0.3.0 official Lightning training path."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from midas_synthetic_training_probe import small_config, synthetic_mudata


def run(output: Path, seed: int) -> dict[str, object]:
    import torch
    import scmidas
    from scmidas import MIDAS

    if output.exists():
        raise RuntimeError("output exists")
    output.mkdir(parents=True, mode=0o750)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    mdata = synthetic_mudata(seed)
    MIDAS.setup_mudata(
        mdata, batch_key="batch", dims_x={"rna": [20], "atac": [10, 20]}
    )
    model = MIDAS(
        mdata,
        configs=small_config(),
        batch_size=16,
        n_save=1000,
        save_model_path=str(output / "unused_checkpoints"),
    )
    model.train(
        max_epochs=1,
        accelerator="gpu",
        devices=1,
        limit_train_batches=2,
        logger=False,
        enable_checkpointing=False,
        enable_progress_bar=False,
        num_sanity_val_steps=0,
    )
    state_path = output / "official_trainer_tensor_state.pt"
    torch.save(model.state_dict(), state_path)
    loaded = torch.load(state_path, map_location="cpu", weights_only=True)
    if set(loaded) != set(model.state_dict()):
        raise RuntimeError("official trainer state census differs")
    receipt = {
        "schema_version": "masld-bench-midas-official-trainer-probe-v1",
        "status": "pass",
        "model_id": "midas",
        "exact_distribution": f"scmidas {scmidas.__version__}",
        "official_trainer_fit_passed": True,
        "epochs": 1,
        "limit_train_batches": 2,
        "safe_weights_only_state_load": True,
        "synthetic_only": True,
        "project_data_read": False,
        "outcomes_read": False,
        "test_outcomes_read": False,
        "champion_claim_allowed": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=8111)
    arguments = parser.parse_args()
    run(arguments.output, arguments.seed)


if __name__ == "__main__":
    main()
