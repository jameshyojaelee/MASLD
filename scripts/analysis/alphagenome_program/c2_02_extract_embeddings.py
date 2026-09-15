#!/usr/bin/env python3
"""C2 step 2: frozen REF/ALT/REF_RC/ALT_RC embeddings for Currin peak-lead SNVs.

Same geometry and pooling as the reporter run: 4,096-bp window with the variant at 0-based index
2048, mean-pooled over positions [1280, 2816).  Model construction, checkpoint restore and the
Caduceus RCPS realignment follow the frozen probes in Analysis/MASLD_Model_Benchmark/scripts
(`extract_gse281364_dna_lm_embeddings.py`, `extract_chromatin_region_dna_lm_embeddings.py`,
`dna_lm_probe_caduceus.py`).

Inputs are token-encoded windows (A7 C8 G9 T10) and the alternate allele token; no measured allele
effect is on disk in anything this script opens.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import random
import sys
import time

import numpy as np
import torch
from safetensors.torch import load_file

BENCH = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/MASLD_Model_Benchmark")
CODE_ROOT = BENCH / "executions/dna-lm-custom-code-21064626/code"
CHECKPOINTS = {
    "hyenadna": BENCH / "executions/dna-lm-open-safetensors-21064556/checkpoints/hyenadna.model.safetensors",
    "caduceus": BENCH / "executions/dna-lm-open-safetensors-21064556/checkpoints/caduceus.model.safetensors",
}
ALLELES = ("REF", "ALT", "REF_RC", "ALT_RC")
WINDOW = 4096
POOL_START0, POOL_END0 = 1280, 2816
WIDTH = 256
SEP_TOKEN = 1


class ExtractionError(RuntimeError):
    """Raised when a model identity, checkpoint restore or pooling requirement differs."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def build_hyenadna(device: torch.device):
    sys.path.insert(0, str(CODE_ROOT))
    from hyenadna.configuration_hyena import HyenaConfig
    from hyenadna.modeling_hyena import HyenaDNAForCausalLM

    code = CODE_ROOT / "hyenadna"
    config = HyenaConfig.from_json_file(str(code / "config.json"))
    if config.d_model != WIDTH or config.max_seq_len != 1_000_002:
        raise ExtractionError("HyenaDNA model identity differs")
    model = HyenaDNAForCausalLM(config)
    state = load_file(str(CHECKPOINTS["hyenadna"]), device="cpu")
    incompatible = model.load_state_dict(state, strict=False)
    aliases = {
        f"hyena.backbone.layers.{layer}.mixer.filter_fn.implicit_filter.{position}.freq"
        for layer in range(config.n_layer)
        for position in (3, 5)
    }
    shared = all(
        model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[1]
        is model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[3]
        is model.hyena.backbone.layers[layer].mixer.filter_fn.implicit_filter[5]
        for layer in range(config.n_layer)
    )
    if set(incompatible.missing_keys) != aliases or incompatible.unexpected_keys or not shared:
        raise ExtractionError("HyenaDNA checkpoint restore differs")
    del state
    model.to(device).eval()

    def embed(tokens: torch.Tensor) -> np.ndarray:
        separator = torch.full((tokens.shape[0], 1), SEP_TOKEN, dtype=torch.long, device=device)
        ids = torch.cat([tokens, separator], dim=1)
        with torch.inference_mode():
            hidden = model.hyena(input_ids=ids, return_dict=True).last_hidden_state
            if hidden.shape[1] != WINDOW + 1 or hidden.shape[2] != WIDTH:
                raise ExtractionError("HyenaDNA hidden shape differs")
            return hidden[:, POOL_START0:POOL_END0].float().mean(dim=1).cpu().numpy()

    return embed, {
        "shared_parameter_aliases": len(aliases),
        "tokens_including_terminal_sep": WINDOW + 1,
        "n_layer": int(config.n_layer),
    }


def build_caduceus(device: torch.device):
    sys.path.insert(0, str(CODE_ROOT))
    from caduceus.configuration_caduceus import CaduceusConfig
    from caduceus.modeling_caduceus import CaduceusForMaskedLM

    code = CODE_ROOT / "caduceus"
    config = CaduceusConfig.from_json_file(str(code / "config.json"))
    if (
        config.d_model != WIDTH
        or config.n_layer != 16
        or not config.rcps
        or not config.bidirectional
        or not config.bidirectional_weight_tie
        or config.bidirectional_strategy != "add"
    ):
        raise ExtractionError("Caduceus model identity differs")
    model = CaduceusForMaskedLM(config)
    state = load_file(str(CHECKPOINTS["caduceus"]), device="cpu")
    incompatible = model.load_state_dict(state, strict=False)
    aliases = {
        f"caduceus.backbone.layers.{layer}.mixer.submodule.mamba_rev.{projection}.weight"
        for layer in range(config.n_layer)
        for projection in ("in_proj", "out_proj")
    } | {"lm_head.lm_head.weight"}
    if set(incompatible.missing_keys) != aliases or incompatible.unexpected_keys:
        raise ExtractionError("Caduceus checkpoint restore differs")
    del state
    tied = all(
        layer.mixer.submodule.mamba_rev.in_proj.weight is layer.mixer.submodule.mamba_fwd.in_proj.weight
        and layer.mixer.submodule.mamba_rev.out_proj.weight is layer.mixer.submodule.mamba_fwd.out_proj.weight
        for layer in model.caduceus.backbone.layers
    ) and model.lm_head.weight is model.get_input_embeddings().weight
    if not tied:
        raise ExtractionError("Caduceus tied parameter identities differ")
    model.to(device).eval()

    def embed(tokens: torch.Tensor) -> np.ndarray:
        with torch.inference_mode():
            hidden = model.caduceus(input_ids=tokens, return_dict=True).last_hidden_state
            if tuple(hidden.shape[1:]) != (WINDOW, 2 * WIDTH):
                raise ExtractionError("Caduceus RCPS hidden shape differs")
            realigned = (hidden[..., :WIDTH] + torch.flip(hidden[..., WIDTH:], dims=(-2, -1))) / 2.0
            return realigned[:, POOL_START0:POOL_END0].float().mean(dim=1).cpu().numpy()

    return embed, {
        "rcps_realignment": "(forward + flip(reverse)) / 2",
        "tokens": WINDOW,
        "n_layer": int(config.n_layer),
    }


BUILDERS = {"hyenadna": build_hyenadna, "caduceus": build_caduceus}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tokens", type=Path, required=True)
    parser.add_argument("--meta", type=Path, required=True)
    parser.add_argument("--model", required=True, choices=sorted(BUILDERS))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260914)
    parser.add_argument("--batch", type=int, default=16)
    parser.add_argument("--limit", type=int, default=0, help="smoke-test row cap; 0 means all rows")
    arguments = parser.parse_args()
    if arguments.output.exists():
        raise ExtractionError(f"refusing to overwrite {arguments.output}")

    random.seed(arguments.seed)
    np.random.seed(arguments.seed)
    torch.manual_seed(arguments.seed)
    torch.cuda.manual_seed_all(arguments.seed)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise ExtractionError("CUDA is required")

    tokens = np.load(arguments.tokens)
    with np.load(arguments.meta, allow_pickle=False) as meta:
        row_index = meta["row_index"]
        variant_id = meta["variant_id"].astype(str)
        ref_token = meta["ref_token"]
        alt_token = meta["alt_token"]
        variant_index0 = meta["variant_index0"]
    n = len(tokens) if arguments.limit <= 0 else min(int(arguments.limit), len(tokens))
    centre = int(variant_index0[0])
    if (
        tokens.shape[1] != WINDOW
        or tokens.dtype != np.uint8
        or len(row_index) != len(tokens)
        or not np.array_equal(tokens[:, centre], ref_token)
        or not np.array_equal(variant_index0, np.full(len(tokens), centre))
        or int(tokens.min()) < 7
        or int(tokens.max()) > 10
    ):
        raise ExtractionError("token window contract differs")

    embed, native = BUILDERS[arguments.model](device)
    embeddings = np.empty((n, len(ALLELES), WIDTH), dtype=np.float32)
    first_value = None
    first_tokens = None
    clock = time.time()
    for begin in range(0, n, arguments.batch):
        stop = min(begin + arguments.batch, n)
        block = torch.from_numpy(tokens[begin:stop].astype(np.int64)).to(device)
        alternative = block.clone()
        alternative[:, centre] = torch.from_numpy(alt_token[begin:stop].astype(np.int64)).to(device)
        # Reverse complement in token space: A7<->T10 and C8<->G9, so token' = 17 - token.
        forms = (block, alternative, 17 - torch.flip(block, dims=(-1,)), 17 - torch.flip(alternative, dims=(-1,)))
        for allele_index, form in enumerate(forms):
            values = embed(form)
            if values.shape != (stop - begin, WIDTH) or not np.isfinite(values).all():
                raise ExtractionError("pooled embedding differs")
            embeddings[begin:stop, allele_index] = values
        if first_value is None:
            first_value = embeddings[begin, 0].copy()
            first_tokens = block[:1].clone()
        if (begin // arguments.batch) % 200 == 0:
            rate = stop / max(1e-9, time.time() - clock)
            print(
                f"{arguments.model} {stop}/{n} variants  {rate:.2f} variants/s  "
                f"eta {(n - stop) / max(rate, 1e-9) / 60:.1f} min",
                flush=True,
            )
    repeated = embed(first_tokens)[0]
    reproducibility = float(np.max(np.abs(repeated - first_value)))
    strand_gap = float(np.max(np.abs(embeddings[:, 0] - embeddings[:, 2])))

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        arguments.output,
        row_index=row_index[:n],
        variant_id=variant_id[:n],
        allele_order=np.asarray(ALLELES),
        embeddings=embeddings,
    )
    receipt = {
        "schema_version": "agp-c2-currin-allele-embeddings-v1",
        "model_id": arguments.model,
        "checkpoint": str(CHECKPOINTS[arguments.model]),
        "checkpoint_sha256": sha256(CHECKPOINTS[arguments.model]),
        "variants": int(n),
        "alleles": list(ALLELES),
        "window_bp": WINDOW,
        "variant_index0": centre,
        "pool_start0": POOL_START0,
        "pool_end0": POOL_END0,
        "pool_length_bp": POOL_END0 - POOL_START0,
        "pooling": "mean_over_pool_window_positions",
        "hidden_width": WIDTH,
        "native_contract": native,
        "reproducibility_max_abs_diff": reproducibility,
        "forward_vs_reverse_complement_max_abs_diff": strand_gap,
        "seed": arguments.seed,
        "batch_size": arguments.batch,
        "device": str(device),
        "torch": torch.__version__,
        "measured_allele_effect_read": False,
        "head_fit": False,
        "elapsed_seconds": round(time.time() - clock, 1),
    }
    Path(str(arguments.output) + ".receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
