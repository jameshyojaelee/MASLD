"""Keep per-library and per-donor cell counts in separate metadata columns."""
import csv
from pathlib import Path

root = Path(__file__).resolve().parents[4]
meta = root / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/metadata"
review = root / "Analysis/SingleCell/reviews/control-label-review-20261021/library_labels.tsv"
with review.open() as f:
    counts = {r["sample"]: int(r["n_cells"]) for r in csv.DictReader(f, delimiter="\t")}
with (meta / "library_metadata.tsv").open() as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
for r in rows:
    r["n_donor_cells"] = r["n_cells"]
    r["n_cells"] = counts[r["sample"]]
assert len(rows) == 275 and sum(r["n_cells"] for r in rows) == 1232318
with (meta / "library_metadata_corrected.tsv").open("x", newline="") as f:
    w = csv.DictWriter(f, rows[0].keys(), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
