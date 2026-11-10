"""Read the source-qualified RNA library labels; never infer disease from dataset."""
import csv
import os
from pathlib import Path


def load_library_diagnoses(base):
    root = Path(os.environ.get("SINGLECELL_DIAGNOSIS_ROOT", str(Path(base) / "Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/metadata")))
    path = root / "library_metadata_corrected.tsv"
    if (root / "library_metadata_source_F.tsv").exists():
        path = root / "library_metadata_source_F.tsv"
    if not path.exists():
        path = root / "library_metadata.tsv"
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {r["sample"]: r for r in rows}
    if len(result) != len(rows):
        raise ValueError("Conflicting RNA library identities")
    return result
