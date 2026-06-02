"""Shared runner base — boilerplate every per-model runner uses.

Per-model runners (state_runner.py, scgpt_runner.py, etc.) instantiate a
ModelAdapter subclass and call `run_arm()`. The base handles:
  - CLI parsing (--modality, --context, --hits, --out)
  - Loading the hit set
  - Loading the substrate (hepatocyte atlas, bulk, or pseudobulk)
  - Iterating hits × contexts and calling `adapter.predict_for_hit(...)`
  - Emitting a validated ModelRunOutput JSON

Each subagent only has to:
  1. Subclass ModelAdapter
  2. Fill `load_checkpoint()` + `predict_for_hit()`  (TODO blocks)
  3. Run `python <model>_runner.py --modality zero_shot --context F__Steatohepatitis__Progressor`
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

# Make `shared/` importable from arm subdirs.
SHARED_DIR = Path(__file__).resolve().parent
if str(SHARED_DIR) not in sys.path:
    sys.path.insert(0, str(SHARED_DIR))

import pandas as pd  # noqa: E402

from output_schema import ModelRunOutput  # noqa: E402

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
PERTURB_ROOT = PROJECT_ROOT / "Analysis/Perturbation"
HITS_DIR = PERTURB_ROOT / "data/hits"
RESULTS_ROOT = PERTURB_ROOT / "results"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

ARM_TO_RESULTS_DIR = {
    "D1": "d1_mechanism",
    "D2": "d2_reversal",
    "D3": "d3_synergy",
    "D4": "d4_circuit",
    "D5": "d5_mouse",
}


def build_arg_parser(model_name: str, arm: str, default_hits: str) -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=f"{model_name} runner for arm {arm}",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--modality",
        choices=["zero_shot", "fine_tuned"],
        required=True,
    )
    p.add_argument(
        "--context",
        default="all",
        help="Stratification context label (e.g. F__Steatohepatitis__Progressor) or 'all'",
    )
    p.add_argument(
        "--hits",
        default=str(HITS_DIR / default_hits),
        help="Path to hit-set CSV",
    )
    p.add_argument(
        "--checkpoint",
        default=None,
        help="Override path to model checkpoint",
    )
    arm_subdir = ARM_TO_RESULTS_DIR.get(arm.upper(), arm.lower())
    p.add_argument(
        "--out-dir",
        default=str(RESULTS_ROOT / arm_subdir),
        help="Output directory for ModelRunOutput JSON",
    )
    p.add_argument(
        "--max-hits",
        type=int,
        default=None,
        help="Cap on hits to score (for debug / dry-run)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Skip model load; emit placeholder predictions only",
    )
    return p


# ---------------------------------------------------------------------------
# Adapter base
# ---------------------------------------------------------------------------


class ModelAdapter:
    """Subclass me. Every per-model runner fills these methods."""

    arm: str = "D?"
    model_name: str = "unknown"
    default_checkpoint: str | None = None

    def __init__(self, modality: str, checkpoint: str | None = None, dry_run: bool = False):
        self.modality = modality
        self.checkpoint = checkpoint or self.default_checkpoint
        self.dry_run = dry_run
        self.model: Any = None
        self.checkpoint_hash: str = "dryrun" if dry_run else ""

    # ---- override --------------------------------------------------------
    def load_checkpoint(self) -> None:
        """TODO: load model weights from self.checkpoint and set self.model."""
        raise NotImplementedError

    def predict_for_hit(
        self,
        *,
        hit_row: pd.Series,
        context: str,
        substrate: Any,
    ) -> list[dict]:
        """TODO: emit a list of *Prediction dicts (matching arm's schema).

        Returned dicts are validated by pydantic when ModelRunOutput is built.
        """
        raise NotImplementedError

    # ---- helpers shared across adapters ---------------------------------
    def _hash_checkpoint(self) -> str:
        if not self.checkpoint:
            return "no_checkpoint"
        p = Path(self.checkpoint)
        if not p.exists():
            return f"missing:{p.name}"
        h = hashlib.sha256()
        if p.is_dir():
            # Hash a deterministic listing of files + sizes (cheap; no I/O on
            # multi-GB safetensors / .ckpt directories). Sufficient for
            # reproducibility tracking — if any file in the checkpoint dir
            # changes size, the hash changes.
            for f in sorted(p.rglob("*")):
                if f.is_file():
                    h.update(f.relative_to(p).as_posix().encode())
                    h.update(b":")
                    try:
                        h.update(str(f.stat().st_size).encode())
                    except OSError:
                        h.update(b"err")
                    h.update(b"\n")
        else:
            with open(p, "rb") as fh:
                while chunk := fh.read(1 << 20):
                    h.update(chunk)
        return h.hexdigest()[:16]


# ---------------------------------------------------------------------------
# Hit loading
# ---------------------------------------------------------------------------


def load_hits(hits_path: str | Path, max_hits: int | None = None) -> pd.DataFrame:
    df = pd.read_csv(hits_path)
    if max_hits is not None:
        df = df.head(max_hits)
    return df


# ---------------------------------------------------------------------------
# Substrate stub — subagents wire real loaders here
# ---------------------------------------------------------------------------


def default_substrate_loader(arm: str, modality: str, context: str) -> Any:
    """Return a (per-arm) substrate handle.

    Concrete loaders live in `shared/load_hepatocyte_atlas.py` +
    `shared/load_bulk_5cohort.py`. This stub returns a small dict
    so a dry-run runner still produces output.
    """
    if arm.upper() == "D1":
        # D1 wants single-cell hepatocyte
        return {"kind": "hepatocyte_atlas", "context": context}
    if arm.upper() == "D2":
        return {"kind": "hepatocyte_atlas+ref_sigs", "context": context}
    if arm.upper() == "D3":
        return {"kind": "hepatocyte_atlas_pairs", "context": context}
    if arm.upper() == "D4":
        return {"kind": "atlas_w_other_celltypes", "context": context}
    if arm.upper() == "D5":
        return {"kind": "mouse_atlas", "context": context}
    return {"kind": "unknown", "context": context}


# ---------------------------------------------------------------------------
# Main run loop
# ---------------------------------------------------------------------------


def run_arm(
    *,
    adapter: ModelAdapter,
    hits_df: pd.DataFrame,
    context: str,
    out_dir: Path,
    substrate_loader=default_substrate_loader,
    notes: str = "",
) -> Path:
    """Drive the prediction loop and emit a validated ModelRunOutput JSON."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    if not adapter.dry_run:
        adapter.load_checkpoint()
        adapter.checkpoint_hash = adapter._hash_checkpoint()

    substrate = substrate_loader(adapter.arm, adapter.modality, context)

    predictions: list[dict] = []
    for _, row in hits_df.iterrows():
        try:
            preds = adapter.predict_for_hit(
                hit_row=row, context=context, substrate=substrate
            )
        except NotImplementedError:
            preds = []
        if adapter.dry_run and not preds:
            preds = _placeholder_predictions(adapter.arm, row)
        predictions.extend(preds)

    elapsed = time.time() - t0
    slurm_job_id = os.environ.get("SLURM_JOB_ID")

    envelope = ModelRunOutput(
        arm=adapter.arm,
        model=adapter.model_name,
        modality=adapter.modality,
        context=context,
        timestamp=datetime.now(timezone.utc),
        checkpoint_hash=adapter.checkpoint_hash or "unknown",
        n_predictions=len(predictions),
        predictions=predictions,
        runtime_seconds=elapsed,
        slurm_job_id=slurm_job_id,
        notes=notes,
    )

    out_path = out_dir / f"{adapter.model_name}_{adapter.modality}_{context}.json"
    with open(out_path, "w") as fh:
        json.dump(envelope.model_dump(mode="json"), fh, indent=2)

    print(f"[{adapter.model_name}/{adapter.modality}] context={context} "
          f"n_pred={len(predictions)} elapsed={elapsed:.1f}s -> {out_path}")
    if slurm_job_id:
        print(f"[{adapter.model_name}] SLURM_JOB_ID={slurm_job_id}")
    return out_path


# ---------------------------------------------------------------------------
# Placeholder predictions so dry-run / skeleton-only mode produces valid JSON.
# ---------------------------------------------------------------------------


def _placeholder_predictions(arm: str, hit_row: pd.Series) -> list[dict]:
    arm = arm.upper()
    if arm == "D1":
        target = str(hit_row.get("gene", "PLACEHOLDER"))
        return [
            dict(
                target_gene=target,
                downstream_gene=f"DOWN_{i}",
                logFC_predicted=0.0,
                abs_rank=i + 1,
                direction="—",
                confidence=0.0,
            )
            for i in range(3)
        ]
    if arm == "D2":
        gene = str(hit_row.get("gene", "PLACEHOLDER"))
        return [
            dict(
                gene=gene,
                reversal_score=0.0,
                reversal_rank=1,
                stage_specific="—",
                cell_type="hepatocyte_progressor",
                reference_signature="consensus",
            )
        ]
    if arm == "D3":
        g1, g2 = (
            str(hit_row.get("gene1", "G1")),
            str(hit_row.get("gene2", "G2")),
        )
        return [
            dict(
                gene1=g1,
                gene2=g2,
                additive_baseline=0.0,
                observed_double_or_higher=0.0,
                synergy_magnitude=0.0,
                synergy_class="additive",
                sigma_above_additive=0.0,
                k562_bias_confidence=0.0,
            )
        ]
    if arm == "D4":
        ligand = str(hit_row.get("ligand", "LIG"))
        rcv_raw = str(hit_row.get("receiver_cell_type", "Macrophage"))
        # Canonicalize to schema's singular Literal options
        _rcv_canon = {
            "macrophages": "Macrophage", "macrophage": "Macrophage",
            "stellates": "Stellate", "stellate": "Stellate", "hsc": "Stellate",
            "lsecs": "LSEC", "lsec": "LSEC", "endothelial": "LSEC",
            "cholangiocytes": "Cholangiocyte", "cholangiocyte": "Cholangiocyte",
        }
        rcv = _rcv_canon.get(rcv_raw.lower(), "Macrophage")
        receptor = str(hit_row.get("receptor", "R"))
        return [
            dict(
                sender_gene=ligand,
                receiver_cell_type=rcv,
                receptor_gene=receptor,
                receiver_response_genes=["R_PLACEHOLDER"],
                response_magnitude=0.0,
                propagation_method="state_2stage",
            )
        ]
    if arm == "D5":
        target = str(hit_row.get("mouse_gene", "Placeholder"))
        return [
            dict(
                target_gene=target,
                downstream_gene=f"down_{i}",
                logFC_predicted=0.0,
                abs_rank=i + 1,
                direction="—",
                confidence=0.0,
            )
            for i in range(3)
        ]
    return []
