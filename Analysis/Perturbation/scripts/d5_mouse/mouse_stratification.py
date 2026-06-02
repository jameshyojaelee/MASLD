"""D5 mouse stratification module.

Generates the (diet × strain × sex) context grid that each D5 runner iterates
over. Pulls from mouse unified metadata.

Subagent: d5-mouse-stratifier.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
MOUSE_META = (
    PROJECT_ROOT / "RNA-seq/Mouse/Unified_Integration/metadata/unified_mouse_metadata.csv"
)


@dataclass
class MouseContext:
    diet_model: str
    strain: str
    sex: str
    n_samples: int


def list_contexts(min_samples: int = 3) -> list[MouseContext]:
    if not MOUSE_META.exists():
        return []
    meta = pd.read_csv(MOUSE_META)
    diet_col = "diet_model" if "diet_model" in meta.columns else "diet"
    sex_col = "sex" if "sex" in meta.columns else "Sex"
    strain_col = "strain" if "strain" in meta.columns else "Strain"
    for c in (diet_col, sex_col, strain_col):
        if c not in meta.columns:
            meta[c] = "unknown"
    grouped = (
        meta.groupby([diet_col, strain_col, sex_col]).size().reset_index(name="n")
    )
    grouped = grouped[grouped["n"] >= min_samples]
    return [
        MouseContext(
            diet_model=str(row[diet_col]),
            strain=str(row[strain_col]),
            sex=str(row[sex_col]),
            n_samples=int(row["n"]),
        )
        for _, row in grouped.iterrows()
    ]


def context_label(ctx: MouseContext) -> str:
    return f"{ctx.diet_model}__{ctx.strain}__{ctx.sex}"
