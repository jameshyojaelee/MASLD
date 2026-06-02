"""D4 circuit runner template.

D4 task: hepatocyte ligand KO -> modified hepatocyte secretome -> downstream
receiver cell-type response (Macrophage / Stellate / LSEC / Cholangiocyte).
Emits CircuitPrediction rows.

Three propagation methods (model_name in the schema doubles as method_id):
  - state_2stage (STATE sender->modified->receiver)
  - niches (NICHES R)
  - commot (COMMOT spatial propagation)

USAGE:
    python state_2stage_runner.py --modality zero_shot --context all
"""
from __future__ import annotations

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

ARM = "D4"
DEFAULT_HITS = "d4_circuit_ligands.csv"


class D4Adapter(ModelAdapter):
    arm = ARM
    model_name = "OVERRIDE_ME"
    propagation_method: str = "state_2stage"
    default_checkpoint: str | None = None

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
        """TODO[d4-runner]: emit CircuitPrediction dicts.

        Schema fields:
          sender_gene, receiver_cell_type, receptor_gene,
          receiver_response_genes (list), response_magnitude, propagation_method.

        Example (state_2stage):
            ligand = hit_row["ligand"]
            receiver = hit_row["receiver_cell_type"]
            modified_hep = self.model.simulate_ko(ligand, hep_substrate)
            response = self.receiver_modelers[receiver].respond(modified_hep)
            return [dict(
                sender_gene=ligand,
                receiver_cell_type=receiver,
                receptor_gene=hit_row.get("receptor", "—"),
                receiver_response_genes=response.top_genes(50),
                response_magnitude=float(response.magnitude),
                propagation_method=self.propagation_method,
            )]
        """
        raise NotImplementedError(
            f"{self.model_name}: implement predict_for_hit()"
        )


def main_for_model(adapter_cls: type[D4Adapter]) -> None:
    parser = build_arg_parser(adapter_cls.model_name, ARM, DEFAULT_HITS)
    args = parser.parse_args()

    hits = load_hits(args.hits, max_hits=args.max_hits)
    adapter = adapter_cls(
        modality=args.modality,
        checkpoint=args.checkpoint,
        dry_run=args.dry_run,
    )
    run_arm(
        adapter=adapter,
        hits_df=hits,
        context=args.context,
        out_dir=Path(args.out_dir),
        notes=f"D4 circuit / propagation={adapter.propagation_method}",
    )
