"""D2 reversal runner — per-model template.

Each per-model runner (scgen_runner.py, cpa_runner.py, state_reversal_runner.py,
tahoe_reversal_runner.py) subclasses `D2Adapter` and customizes:
  - MODEL_NAME, DEFAULT_CHECKPOINT
  - load_checkpoint(self)
  - predict_reversal(self, gene, context, substrate, reference_signature)

D2 task: rank ALL atlas genes (33,943) by reversal score (Diseased→Healthy)
within a given hepatocyte context, against each of 3 reference signatures.

USAGE:
    python scgen_runner.py --modality zero_shot --context F__Steatohepatitis__Progressor \\
        --reference ref_a
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import pandas as pd  # noqa: E402

from runner_base import (  # noqa: E402
    ModelAdapter,
    RESULTS_ROOT,
    build_arg_parser,
    load_hits,
    run_arm,
)

ARM = "D2"
DEFAULT_HITS = "d2_reversal_hits.csv"

REFERENCE_SIGNATURES = ("ref_a", "ref_b", "ref_c")
REF_SIG_DIR = (
    Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
    / "Analysis/Perturbation/data/reference_signatures"
)

# Canonical filenames produced by task #18 (d2-ref-curator).
REF_FILENAMES = {
    "ref_a": "ref_a_subtype.csv",      # hepatocyte subtype Progressor+Moderate vs Healthy+Stable
    "ref_b": "ref_b_bulk5cohort.csv",  # bulk 5-cohort dream Disease-vs-Control
    "ref_c": "ref_c_F0vF4.csv",        # hepatocyte F0 vs F4 (Cirrhosis vs Healthy)
}


def load_reference_signature(ref_label: str) -> pd.DataFrame:
    """Load one of three alternative Diseased→Healthy reference signatures.

    Schema (all three CSVs):
        gene, healthy_mean, diseased_mean, LFC, pval, padj
    NOTE: column is `LFC`, NOT `logFC`. ref_b has NaN for healthy_mean/diseased_mean
    (bulk dream output only ships LFC + AveExpr); use LFC as the projection axis.
    """
    filename = REF_FILENAMES.get(ref_label, f"{ref_label}.csv")
    path = REF_SIG_DIR / filename
    if not path.exists():
        # Fall back to an empty signature with the canonical schema; runner
        # still produces valid JSON (predictions will score zero).
        return pd.DataFrame(columns=["gene", "healthy_mean", "diseased_mean", "LFC", "pval", "padj"])
    return pd.read_csv(path)


class D2Adapter(ModelAdapter):
    arm = ARM
    model_name = "OVERRIDE_ME"
    default_checkpoint: str | None = None

    def __init__(self, *args, reference_signature: str = "ref_a", **kwargs):
        super().__init__(*args, **kwargs)
        self.reference_signature = reference_signature
        self._ref_df = load_reference_signature(reference_signature)

    # ----- TODO: subagent fills these in ---------------------------------
    def load_checkpoint(self) -> None:
        raise NotImplementedError(
            f"{self.model_name}: implement load_checkpoint()"
        )

    def predict_for_hit(
        self,
        *,
        hit_row: pd.Series,
        context: str,
        substrate,
    ) -> list[dict]:
        """TODO[d2-runner]: compute reversal score for hit_row.gene vs reference.

        Schema: ReversalPrediction
          gene, reversal_score, reversal_rank, stage_specific, cell_type, reference_signature.

        Example:
            gene = hit_row["gene"]
            ko_state = self.model.simulate_ko(gene, baseline=substrate)
            # NB: ref CSV column is `LFC` (not `logFC`); see load_reference_signature
            score = -cosine(ko_state - baseline, self._ref_df.set_index("gene")["LFC"])
            return [dict(
                gene=gene,
                reversal_score=score,
                reversal_rank=0,  # rank assigned in post-process below
                stage_specific="—",
                cell_type="hepatocyte_progressor",
                reference_signature=self.reference_signature,
            )]
        """
        raise NotImplementedError(
            f"{self.model_name}: implement predict_for_hit()"
        )


def main_for_model(adapter_cls: type[D2Adapter]) -> None:
    parser = build_arg_parser(adapter_cls.model_name, ARM, DEFAULT_HITS)
    parser.add_argument(
        "--reference",
        choices=REFERENCE_SIGNATURES,
        default="ref_a",
        help="Which Diseased->Healthy reference signature to use",
    )
    args = parser.parse_args()

    hits = load_hits(args.hits, max_hits=args.max_hits)
    print(f"[{adapter_cls.model_name}] {len(hits)} hits, ref={args.reference}")

    adapter = adapter_cls(
        modality=args.modality,
        checkpoint=args.checkpoint,
        dry_run=args.dry_run,
        reference_signature=args.reference,
    )

    out_dir = Path(args.out_dir)
    context_label = f"{args.context}__{args.reference}"
    out_path = run_arm(
        adapter=adapter,
        hits_df=hits,
        context=context_label,
        out_dir=out_dir,
        notes=f"D2 reversal / ref={args.reference}",
    )

    # Post-process: assign reversal_rank from scores so consensus aggregator
    # can rely on the field being populated even when subagents forget.
    import json

    with open(out_path) as fh:
        envelope = json.load(fh)
    preds = envelope["predictions"]
    if preds:
        preds.sort(key=lambda p: p.get("reversal_score", 0.0), reverse=True)
        for i, p in enumerate(preds):
            p["reversal_rank"] = i + 1
        envelope["predictions"] = preds
        with open(out_path, "w") as fh:
            json.dump(envelope, fh, indent=2)
