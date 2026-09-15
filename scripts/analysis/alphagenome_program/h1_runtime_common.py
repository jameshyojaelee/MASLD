#!/usr/bin/env python3
"""Shared loading and scoring helpers for the local AlphaGenome runtime probe.

The local weights are the gated `google/alphagenome-all-folds` Orbax checkpoint acquired at
`executions/alphagenome-weights-20260915T103442Z`. The code is `google-deepmind/alphagenome_research`
pinned at 1e55dcffb98ba26b31e74edc5e9f038f54c0e89d.

The scoring functions here are ported from the hosted-API recipe in
`scripts/analysis/alphagenome_atlas/63_indel_rescue.py` so that the local value and the hosted value are
the same quantity computed the same way. Nothing in that directory is modified or imported.

Outputs of these LOCAL weights may be used as features or teacher signals. Hosted Atlas/API outputs may
NOT be used for training, and the two are never averaged or pooled.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import pathlib
import time

import numpy as np

PROJ = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CHECKPOINT = (
    PROJ
    / "Analysis/MASLD_Model_Benchmark/executions/alphagenome-weights-20260915T103442Z/checkpoints"
)
FASTA_PATH = (
    "/gpfs/commons/home/jameslee/reference_genome/cellranger-atac/"
    "refdata-cellranger-arc-GRCh38-2024-A/fasta/genome.fa"
)

# The hosted recipe's configuration, reproduced exactly (63_indel_rescue.py main()).
OUTPUT_NAMES = ["RNA_SEQ", "ATAC", "DNASE", "CHIP_HISTONE", "SPLICE_SITE_USAGE"]
LIVER_TERMS = ["UBERON:0002107", "UBERON:0001114", "UBERON:0001115", "CL:0000182"]
FLANK = 2_000
EPS = 1e-9

CHANNELS = (
    ("rna", "rna_seq"),
    ("atac", "atac"),
    ("dnase", "dnase"),
    ("h3k27ac", "chip_histone"),
    ("splice", "splice_site_usage"),
)


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}] {msg}", flush=True)


def sha256_file(path: str | os.PathLike[str]) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------------------------------
# model
# --------------------------------------------------------------------------------------------------
def load_model(checkpoint_path=CHECKPOINT, fasta_path: str | None = FASTA_PATH):
    """Load the local checkpoint onto the first GPU.

    `default_organism_settings()` points fasta/gtf/splice-site/calibration at remote Google Storage URLs.
    Those are replaced with the project's own GRCh38 FASTA (human) and with None (everything else), so
    nothing is fetched at load time and no annotation source other than the audited one enters. Metadata
    is left at the packaged default for BOTH organisms, which keeps Orbax's strict checkpoint validation
    on: dropping mouse would change the parameter tree the restore is validated against.
    """
    import dataclasses

    import jax
    from alphagenome.models import dna_model as dna_model_api
    from alphagenome_research.model import dna_model

    settings = dict(dna_model.default_organism_settings())
    human = dna_model_api.Organism.HOMO_SAPIENS
    mouse = dna_model_api.Organism.MUS_MUSCULUS
    settings[human] = dataclasses.replace(
        settings[human],
        fasta_path=fasta_path,
        gtf_feather_path=None,
        pas_feather_path=None,
        splice_site_starts_feather_path=None,
        splice_site_ends_feather_path=None,
        calibration_path=None,
    )
    settings[mouse] = dataclasses.replace(
        settings[mouse],
        fasta_path=None,
        gtf_feather_path=None,
        pas_feather_path=None,
        splice_site_starts_feather_path=None,
        splice_site_ends_feather_path=None,
        calibration_path=None,
    )

    device = jax.local_devices()[0]
    if device.platform != "gpu":
        raise RuntimeError(f"expected a GPU device, jax offers {jax.local_devices()}")
    t0 = time.time()
    model = dna_model.create(
        str(checkpoint_path), organism_settings=settings, device=device
    )
    load_s = time.time() - t0
    log(f"checkpoint restored in {load_s:.1f}s onto {device}")
    return model, device, load_s


def requested_outputs():
    from alphagenome.models import dna_client

    return [getattr(dna_client.OutputType, o) for o in OUTPUT_NAMES]


def memory_stats(device) -> dict:
    try:
        st = device.memory_stats() or {}
    except Exception as exc:  # pragma: no cover - platform dependent
        return {"error": repr(exc)}
    return {k: int(v) for k, v in st.items() if isinstance(v, (int, float))}


# --------------------------------------------------------------------------------------------------
# sequence construction (ported verbatim in behaviour from 63_indel_rescue.py)
# --------------------------------------------------------------------------------------------------
class ContractError(RuntimeError):
    pass


def indel_alt_sequence(fetch, chrom: str, start0: int, end: int, pos1: int, ref: str, alt: str) -> str:
    """Alternate window, same length as the reference window, with the edit applied at pos1.

    An insertion pushes bases off the far end; a deletion pulls bases in from beyond it. The readout is
    local to the variant at the window centre, so the compensated edge is ~500 kb away.
    """
    width = end - start0
    i = pos1 - 1 - start0
    if i < 0 or i + len(ref) > width:
        raise ContractError(f"variant {chrom}:{pos1} {ref}>{alt} does not fit inside [{start0},{end})")
    extended = fetch(chrom, start0, end + max(0, len(ref) - len(alt))).upper()
    observed = extended[i : i + len(ref)]
    if observed != ref.upper():
        raise ContractError(
            f"reference mismatch at {chrom}:{pos1}: window has {observed}, variant says {ref}"
        )
    out = extended[:i] + alt.upper() + extended[i + len(ref) :]
    if len(out) < width:
        raise ContractError(f"cannot compensate window length for {chrom}:{pos1}")
    return out[:width]


def fold_mask(mask, n_rows: int):
    """Fold a base-resolution mask onto a track's own row count (CHIP_HISTONE is 128-bp binned)."""
    m = np.asarray(mask, dtype=bool)
    if m.size == n_rows:
        return m
    if n_rows <= 0 or m.size % n_rows:
        raise ContractError(f"mask of {m.size} does not fold onto {n_rows} rows")
    fold = m.size // n_rows
    return m.reshape(n_rows, fold).any(axis=1)


def h3k27ac_columns(track_names) -> list:
    return [("H3K27AC" in str(n).upper()) for n in track_names]


def summarise(output, mask_local) -> dict:
    """Mean |track| over the masked rows, per channel. Identical to the hosted recipe's summarise()."""
    out = {}
    for name, attr in CHANNELS:
        td = getattr(output, attr, None)
        if td is None or td.values is None or td.values.size == 0:
            out[name] = math.nan
            continue
        v = np.asarray(td.values, dtype=float)
        keep = np.ones(v.shape[1], bool)
        if name == "h3k27ac":
            md = getattr(td, "metadata", None)
            if md is None or "name" not in md:
                out[name] = math.nan
                continue
            keep = np.asarray(h3k27ac_columns(md["name"].astype(str)), bool)
            if not keep.any():
                out[name] = math.nan
                continue
        m = fold_mask(mask_local, v.shape[0])
        out[name] = float(np.nanmean(np.abs(v[m, :][:, keep]))) if m.any() else math.nan
    return out


def log2_ratio(alt: float, ref: float) -> float:
    if alt != alt or ref != ref:
        return math.nan
    return math.log2((alt + EPS) / (ref + EPS))


def center_mask(width: int, pos1: int, start0: int, flank: int = FLANK):
    mask = np.zeros(width, bool)
    lo, hi = max(0, pos1 - 1 - start0 - flank), min(width, pos1 - 1 - start0 + flank)
    mask[lo:hi] = True
    return mask


def write_json(path, obj) -> None:
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        json.dump(obj, fh, indent=1, default=_default)
    log(f"wrote {path}")


def _default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)
