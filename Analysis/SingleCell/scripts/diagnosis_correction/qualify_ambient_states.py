"""Condition attenuation labels on the corrected primary GLS family."""
import sys
from pathlib import Path

import pandas as pd

root = Path(__file__).resolve().parents[4]
task = root / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
sys.path.insert(0, str(root / "Analysis/SingleCell/scripts/ambient_program_recalibration"))
from common import evidence_state

effects = pd.read_csv(task / "OLS_ambient/ambient_program_effects.tsv", sep="\t")
primary = pd.read_csv(task / "programs/results_primary_family.tsv", sep="\t").set_index(["cell_type", "module"])
assert len(primary) == 117
for i, r in effects.iterrows():
    p = primary.loc[(r.cell_type, r.module)]
    effects.loc[i, "primary_selected_frozen"] = bool(p.selected_H59)
    effects.loc[i, "robust_display_frozen"] = bool(p.robust_display_H59)
effects["primary_selection_model"] = "Corrected GLS; within-dataset 19999 permutations; BH over 117"
common = effects.analysis_universe.eq("complete_case_common_universe")
effects.loc[common, "evidence_state"] = effects.loc[common].apply(evidence_state, axis=1)
effects.loc[~common, "evidence_state"] = "sensitivity_only_includes_unresolved_GSE189600_records"
out = task / "ambient_program_effects.tsv"
assert not out.exists()
effects.to_csv(out, sep="\t", index=False)
