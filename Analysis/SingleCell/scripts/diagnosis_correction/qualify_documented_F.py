"""Separate deposited numeric SAF from diagnoses and free-text control entries."""
import csv
import re
from pathlib import Path

root = Path(__file__).resolve().parents[4]
task = root / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026"
meta = task / "metadata"


def read(p):
    with p.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write(p, rows):
    with p.open("x", newline="") as f:
        w = csv.DictWriter(f, rows[0].keys(), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


source = read(root / "Analysis/SingleCell/reviews/control-label-review-20261021/primary_source_diagnosis_checks.tsv")
by_donor = {}
for r in source:
    if r["dataset"] != "GSE202379":
        continue
    match = re.fullmatch(r"S\d+A\d+F([0-4])", r["deposited_histology"])
    value = match.group(1) if match else ""
    if r["donor"] in by_donor:
        assert by_donor[r["donor"]]["F_stage_documented"] == value
    by_donor[r["donor"]] = dict(donor=r["donor"], dataset=r["dataset"], F_stage_documented=value,
        deposited_SAF=r["deposited_histology"], source_accession=r["gsm"])
assert by_donor["GSE202379_P98"]["F_stage_documented"] == "1"
assert all(by_donor[d]["F_stage_documented"] == "" for d in ("GSE202379_PHL1", "GSE202379_PHL2"))
write(meta / "documented_F_source.tsv", list(by_donor.values()))
rows = read(meta / "library_metadata_corrected.tsv")
for r in rows:
    r["F_stage_documented"] = by_donor.get(r["donor"], {}).get("F_stage_documented", "")
write(meta / "library_metadata_source_F.tsv", rows)
rows = read(meta / "donor_metadata_stage_corrected.tsv")
for r in rows:
    if r["dataset"] == "GSE202379":
        r["F_stage_documented"] = by_donor[r["donor"]]["F_stage_documented"]
        r["F_stage_doc_source"] = "deposited numeric SAF only; free-text healthy control is not F0"
write(meta / "donor_metadata_stage_corrected_source_F.tsv", rows)
for tag in ("raw", "corrected"):
    rows = read(task / "N3" / tag / "donor_meta.tsv")
    for r in rows:
        if r["dataset"] == "GSE202379":
            r["F_documented"] = by_donor[r["donor"]]["F_stage_documented"]
    write(task / "N3" / tag / "donor_meta_source_F.tsv", rows)
print("Numeric SAF donors:", sum(bool(r["F_stage_documented"]) for r in by_donor.values()))
