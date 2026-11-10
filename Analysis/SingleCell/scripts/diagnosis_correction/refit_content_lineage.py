"""Refit the original content-lineage tests on unchanged cached scores."""
import json
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[4]
task = root / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
source = root / "scripts/analysis/fig34_corrections/03_lineage.py"
sys.path.insert(0, str(source.parent))
import common

old = root / "figures/candidates/fig34-scientific-corrections-20260929-v1/lineage_matched"
scores = common.read(old / "lineage_scores.tsv.gz")
disjoint = common.read(old / "disjoint_lineage_scores.tsv.gz")
text = source.read_text()
body = text[text.index("    rows=[]\n"):text.index("    snapshot([")]
body = textwrap.dedent(body)
target = {list(common.FOCAL)[0]: "Fibroblasts", list(common.FOCAL)[1]: "Cholangiocytes"}


def calculate(s, dis, out):
    out.mkdir(parents=True, exist_ok=False)
    namespace = dict(vars(common), s=s, dis=dis, out=out, TARGET=target)
    exec(compile(body, str(source), "exec"), namespace)
    return namespace["fits"]


reproduced = calculate(scores, disjoint, task / "content_lineage_original_reproduction")
original = common.read(old / "lineage_cohort_effects.tsv")
keys = ["program_uid", "dataset", "minimum_cells", "annotation_filter", "arm", "contrast"]
a = original.sort_values(keys).reset_index(drop=True)
b = reproduced.sort_values(keys).reset_index(drop=True)
assert a[keys].equals(b[keys])
assert np.allclose(a[["beta", "hc3_se", "hc3_p", "hc3_q"]],
                   b[["beta", "hc3_se", "hc3_p", "hc3_q"]], equal_nan=True, atol=1e-10)
metadata = common.read(task / "metadata/donor_metadata.tsv").set_index("donor")
out = task / "content_lineage"
for frame in [scores, disjoint]:
    frame["source_diagnosis_previous"] = frame.disease
    frame["disease"] = frame.donor.map(metadata.disease_stage_coarse)
    frame["source_excluded"] = frame.donor.map(metadata.exclude_stage_analysis)
    assert frame.disease.notna().all() and frame.source_excluded.notna().all()
calculate(scores[~scores.source_excluded].copy(), disjoint[~disjoint.source_excluded].copy(), out)
common.write(scores, out / "lineage_scores_with_eligibility.tsv.gz")
common.write(disjoint, out / "disjoint_lineage_scores_with_eligibility.tsv.gz")
common.snapshot([source, old / "lineage_scores.tsv.gz", old / "disjoint_lineage_scores.tsv.gz",
                 task / "metadata/donor_metadata.tsv"], out / "inputs.tsv")
(out / "definition.txt").write_text(
    "Fixed projected scores and raw gene calibration unchanged. Original cohort fits reproduced.\n"
    "Corrected deposited RNA diagnoses; P98 excluded. Donor-level MASH-minus-steatosis contrasts.\n"
    "HC3 cohort estimates; native inverse-variance meta-analysis; BH over two programs in each declared arm.\n"
    "50-cell primary; 20-cell, confidence and full-raw sensitivity arms retained. Seed 20260929.\n"
    "Captured-lineage fractions are not whole-tissue cell proportions.\n")
(out / "environment.json").write_text(json.dumps(dict(python=sys.version, numpy=np.__version__, pandas=pd.__version__), indent=2))
(task / "CONTENT_LINEAGE_COMPLETED").touch()
