#!/usr/bin/env python
"""Direct ST-head training launcher.

Bypasses the broken `state tx train` Hydra CLI (arc-state 0.10.2 ships
`defaults: - wandb: default` in its packaged config.yaml but the
`configs/wandb/default.yaml` file itself is missing → MissingConfigException
before any of our overrides apply).

This launcher composes the same DictConfig that `state tx train` would have
composed, fills in the missing `wandb` section locally, and dispatches into
`state._cli._tx._train.run_tx_train(cfg)` directly. Same code path otherwise,
including pl.seed_everything, datamodule construction, and trainer.fit.

Run from inside `perturbation_state` with
  PYTHONPATH=...vendor PYTHONNOUSERSITE=1 python train_st_head_v1.py
(see train_st_head_v1.sbatch).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import yaml
from omegaconf import OmegaConf

# Import after env vars are set externally (sbatch wraps this).
from state._cli._tx._train import run_tx_train


def build_cfg(args: argparse.Namespace) -> dict:
    """Compose the run cfg by-hand (mirrors arc-state 0.10.2 config.yaml shape).

    Keep the same key names as the published `state-defaults.yaml` and the
    `st-se-replogle-full/hepg2_0.99/config.yaml` reference so downstream
    `state tx predict` can re-read the resulting run dir without surprises.
    """
    cfg = {
        "data": {
            "name": "PerturbationDataModule",
            "kwargs": {
                "toml_config_path": args.toml,
                "embed_key": args.embed_key,
                "output_space": args.output_space,
                "pert_rep": "onehot",
                "basal_rep": "sample",
                "num_workers": args.num_workers,
                "pin_memory": True,
                "n_basal_samples": 1,
                "basal_mapping_strategy": "random",
                "should_yield_control_cells": True,
                "batch_col": "gem_group",
                "pert_col": "gene",
                "cell_type_key": "cell_type",
                "control_pert": "non-targeting",
                "map_controls": True,
                "perturbation_features_file": None,
                "store_raw_basal": False,
                "int_counts": False,
                "barcode": True,
            },
            "output_dir": None,
            "debug": False,
        },
        "model": {
            "name": "state",
            "checkpoint": None,
            "device": "cuda",
            "kwargs": {
                "cell_set_len": 64,
                "blur": 0.05,
                # NOTE: published hepg2_0.99 used hidden_dim=328 with
                # num_attention_heads=12 — modern transformers (with the
                # huggingface_hub StrictDataclass validator) reject 328%12 != 0
                # at LlamaConfig construction time. Bump to 768 (12 * 64 head_dim)
                # to match the state.yaml default and stay within strict validators.
                "hidden_dim": 768,
                "loss": "energy",
                "confidence_token": False,
                "n_encoder_layers": 1,
                "n_decoder_layers": 1,
                "predict_residual": True,
                "softplus": True,
                "freeze_pert_backbone": False,
                "transformer_decoder": False,
                "finetune_vci_decoder": False,
                "residual_decoder": False,
                "batch_encoder": True,
                "use_batch_token": False,
                "nb_decoder": False,
                "mask_attn": False,
                "use_effect_gating_token": False,
                "distributional_loss": "energy",
                "init_from": None,
                "mmd_num_chunks": 1,
                "randomize_mmd_chunks": False,
                "llm_name": None,
                "transformer_backbone_key": "llama",
                "transformer_backbone_kwargs": {
                    "bidirectional_attention": True,
                    "max_position_embeddings": 64,
                    "hidden_size": 768,
                    "intermediate_size": 3072,
                    "num_hidden_layers": 8,
                    "num_attention_heads": 12,
                    "num_key_value_heads": 12,
                    "head_dim": 64,
                    "use_cache": False,
                    "attention_dropout": 0.0,
                    "hidden_dropout": 0.0,
                    "layer_norm_eps": 1.0e-06,
                    "pad_token_id": 0,
                    "bos_token_id": 1,
                    "eos_token_id": 2,
                    "tie_word_embeddings": False,
                    "rotary_dim": 0,
                    "use_rotary_embeddings": False,
                },
                "lora": {
                    "enable": False,
                    "r": 16,
                    "alpha": 32,
                    "dropout": 0.05,
                    "bias": "none",
                    "target": "auto",
                    "adapt_mlp": False,
                    "task_type": "FEATURE_EXTRACTION",
                    "merge_on_eval": False,
                },
            },
        },
        "training": {
            "wandb_track": False,
            "weight_decay": 0.0005,
            "batch_size": args.batch_size,
            "lr": 0.0001,
            "max_steps": args.max_steps,
            "train_seed": 42,
            "val_freq": args.val_freq,
            "ckpt_every_n_steps": args.ckpt_every_n_steps,
            "gradient_clip_val": 10,
            "gradient_accumulation_steps": 1,
            "loss_fn": "mse",
            "devices": 1,
            "strategy": "auto",
            "use_mfu": True,
            "mfu_kwargs": {
                "available_flops": 60_000_000_000_000.0,
                "use_backward": True,
                "logging_interval": 10,
                "window_size": 2,
            },
            "cumulative_flops_use_backward": True,
        },
        # Required by run_tx_train even when use_wandb=False (it reads
        # cfg.wandb.{entity,project,local_wandb_dir} into get_loggers).
        "wandb": {
            "entity": "nslab",
            "project": "state-masld",
            "local_wandb_dir": str(Path(args.output_dir).parent / "wandb_logs"),
            "tags": ["masld_st_head"],
        },
        "name": args.name,
        "output_dir": args.output_dir,
        "use_wandb": False,
        "overwrite": args.overwrite,
        "return_adatas": False,
        "pred_adata_path": None,
        "true_adata_path": None,
    }
    return cfg


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--toml", required=True, help="cell-load TOML config path")
    ap.add_argument("--name", required=True, help="run name (subdir of --output-dir)")
    ap.add_argument("--output-dir", required=True, help="checkpoint root")
    ap.add_argument("--embed-key", default="X_hvg")
    ap.add_argument("--output-space", default="gene")
    ap.add_argument("--num-workers", type=int, default=8)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--max-steps", type=int, default=80_000)
    ap.add_argument("--val-freq", type=int, default=2_000)
    ap.add_argument("--ckpt-every-n-steps", type=int, default=4_000)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    cfg_dict = build_cfg(args)
    print("[train-st-v1] cfg yaml head:")
    print(yaml.safe_dump({k: cfg_dict[k] for k in ("name", "output_dir", "use_wandb", "wandb")}, default_flow_style=False))

    # run_tx_train expects a DictConfig (it calls OmegaConf.to_container).
    cfg = OmegaConf.create(cfg_dict)
    run_tx_train(cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
