"""Read-only review of atlas diagnosis labels; never opens expression matrices.

Run with the existing rapids_singlecell Python environment. Outputs are review
evidence only; source labels, results, and manuscript authorities are untouched.
"""
import collections
import csv
import hashlib
import json
import platform
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import openpyxl

ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parent
COR = ROOT / "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014"
inputs = []


def read_table(path, delimiter="\t"):
    inputs.append(path)
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write_table(name, rows):
    if not rows:
        return
    with (OUT / name).open("x", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


run_donor = {}
run_source = {}
for ds in ("GSE244832", "GSE202379", "GSE185477", "GSE136103"):
    for row in read_table(ROOT / f"data/{ds}/metadata/donor_pairing.csv", ","):
        for run in row["rna_srrs"].split(";"):
            run = run.strip()
            assert run not in run_donor
            run_donor[run] = ds + "_" + row["donor_id"]
            run_source[run] = row["condition"]

la_samples = read_table(ROOT / "Liver_Atlas/metadata/GSE192740_sampleInfo_scRNAseq.tsv")
la_short = {r["characteristics: shortFileName"]: r for r in la_samples}
for row in read_table(ROOT / "Liver_Atlas/metadata/SraRunTable.csv", ","):
    sample = la_short.get(row["shortfilename"])
    if sample:
        donor = re.search(r"H\d+", sample["title"])
        if donor:
            assert row["Run"] not in run_donor
            run_donor[row["Run"]] = "Liver_Atlas_" + donor[0]

roster = {r["donor"]: r for r in read_table(COR / "metadata/current_source_diagnosis_donor_roster.tsv")}
source_join = {r["donor"]: r for r in read_table(COR / "source/Liver_Atlas_TableS6_donor_join.tsv")}
xlsx = COR / "source/Guilliams_TableS6.xlsx"
inputs.append(xlsx)
wb = openpyxl.load_workbook(xlsx, read_only=True, data_only=True)
sheet = wb["Human Samples"]
records = list(sheet.values)
header = records[0][:13]
histology = {r[0]: dict(zip(header, r[:13])) for r in records[1:] if r[0]}
wb.close()
for donor, row in source_join.items():
    native = histology[row["source_sample"]]
    assert str(native["Steatosis"]) == row["source_steatosis"]
    assert str(native["Fibrosis"]) == row["source_fibrosis"]

# Only categorical obs codes and their category names are read, about 10 MB.
# No X, raw/X, layers, embeddings, graphs, or gene counts are opened.
atlas = ROOT / "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
columns = ["sample", "dataset", "condition", "condition_harmonized", "preparation_method"]
with h5py.File(atlas, "r") as handle:
    categories = [handle[f"obs/{c}/categories"].asstr()[:] for c in columns]
    arrays = [handle[f"obs/{c}/codes"][:] for c in columns]
    assert all(np.all(a >= 0) for a in arrays)
    keys, counts = np.unique(np.stack(arrays, axis=1), axis=0, return_counts=True)
    cell_types = handle["obs/cell_type/categories"].asstr()[:]
    ct_codes = handle["obs/cell_type/codes"][:]
    assert np.all(ct_codes >= 0)
    ct_counts = np.bincount(
        arrays[0].astype(np.int64) * len(cell_types) + ct_codes,
        minlength=len(categories[0]) * len(cell_types),
    ).reshape(len(categories[0]), len(cell_types))

assert len(keys) == len(categories[0]), "Within-library metadata conflict"
samples = []
for codes, count in zip(keys, counts):
    row = {c: str(categories[i][codes[i]]) for i, c in enumerate(columns)}
    run = row["sample"]
    donor = run_donor.get(run, run)
    source = roster.get(donor, {})
    stage = source.get("disease_stage_coarse", "Unresolved")
    source_binary = "Healthy" if stage == "Healthy" else (
        "Disease" if stage in ("Steatosis", "Steatohepatitis", "Cirrhosis") else "Unresolved"
    )
    old_binary = "Disease" if row["condition_harmonized"] == "MASLD" else row["condition_harmonized"]
    row.update(
        donor=donor, n_cells=int(count), source_stage=stage,
        source_binary=source_binary, pairing_diagnosis=run_source.get(run, ""),
        source_stage_excluded=source.get("exclude_stage_analysis", "not_in_102_donor_roster"),
        analyzed_donor=bool(source),
        binary_label_problem=(old_binary != source_binary),
        prep_eligible=row["preparation_method"] in ("nuclei", "cd45_negative", "unsorted"),
    )
    for j, name in enumerate(cell_types):
        row["n_" + name.replace(" ", "_")] = int(ct_counts[codes[0], j])
    samples.append(row)
write_table("library_labels.tsv", samples)

donors = []
for donor in sorted({r["donor"] for r in samples}):
    runs = [r for r in samples if r["donor"] == donor]
    for col in ("dataset", "condition", "condition_harmonized", "source_stage"):
        assert len({r[col] for r in runs}) == 1, (donor, col)
    first = runs[0]
    row = {k: first[k] for k in ("donor", "dataset", "condition", "condition_harmonized", "source_stage", "source_binary", "source_stage_excluded", "analyzed_donor", "binary_label_problem")}
    eligible = [r for r in runs if r["prep_eligible"]]
    row.update(n_libraries=len(runs), n_cells=sum(r["n_cells"] for r in runs),
               n_prep_eligible_libraries=len(eligible),
               libraries=";".join(r["sample"] for r in runs),
               prep_eligible_libraries=";".join(r["sample"] for r in eligible))
    native = histology.get(donor.removeprefix("Liver_Atlas_"), {}) if first["dataset"] == "Liver_Atlas" else {}
    for field in ("Steatosis", "Fibrosis", "Reason for hospitalization", "Procedure", "Alcohol intake"):
        row["native_" + field.replace(" ", "_")] = native.get(field, "")
    for name in cell_types:
        col = "n_" + name.replace(" ", "_")
        row[col] = sum(r[col] for r in runs)
    donors.append(row)
write_table("donor_labels.tsv", donors)
write_table("analyzed_donor_label_problems.tsv", [r for r in donors if r["analyzed_donor"] and r["binary_label_problem"]])

summaries = []
for ds in sorted({r["dataset"] for r in donors}):
    rows = [r for r in donors if r["dataset"] == ds]
    summaries.append(dict(dataset=ds, libraries=sum(r["n_libraries"] for r in rows),
                          cells=sum(r["n_cells"] for r in rows), donor_identifiers=len(rows),
                          analyzed_donor_identifiers=sum(r["analyzed_donor"] for r in rows),
                          analyzed_binary_label_problems=sum(r["binary_label_problem"] and r["analyzed_donor"] for r in rows),
                          atlas_healthy_donors=sum(r["condition_harmonized"] == "Healthy" for r in rows),
                          source_healthy_donors=sum(r["source_stage"] == "Healthy" for r in rows)))
write_table("dataset_summary.tsv", summaries)

# Reconstruct the labels eligible for the canonical binary DE script, without
# loading a count matrix or refitting a model. Presence of a pseudobulk column
# is checked from its CSV header only, as in the existing matched-donor step.
pb = ROOT / "Analysis/SingleCell/results_gpu_v2/pseudobulk"
eligible = {r["donor"]: r for r in donors if r["n_prep_eligible_libraries"] and r["n_cells"] >= 50}
membership = []
for path in sorted(pb.glob("*_pseudobulk.csv")):
    if " " in path.name:
        continue
    inputs.append(path)
    with path.open() as handle:
        names = next(csv.reader(handle))[1:]
    present = {run_donor.get(s, s) for s in names}
    for donor in sorted(present & eligible.keys()):
        row = eligible[donor]
        membership.append(dict(cell_type=path.name.removesuffix("_pseudobulk.csv"),
                               donor=donor, dataset=row["dataset"],
                               atlas_binary=row["condition_harmonized"],
                               source_stage=row["source_stage"], source_binary=row["source_binary"],
                               binary_label_problem=row["binary_label_problem"]))
write_table("binary_DE_reconstructed_membership.tsv", membership)

summary = dict(
    atlas_path=str(atlas.resolve()), atlas_obs_cells=int(sum(counts)),
    integrated_libraries=len(samples), integrated_donor_identifiers=len(donors),
    analyzed_donor_identifiers=sum(r["analyzed_donor"] for r in donors),
    analyzed_libraries=sum(r["n_libraries"] for r in donors if r["analyzed_donor"]),
    analyzed_cells=sum(r["n_cells"] for r in donors if r["analyzed_donor"]),
    within_library_label_conflicts=0, within_donor_label_conflicts=0,
    analyzed_atlas_healthy_donors=sum(r["analyzed_donor"] and r["condition_harmonized"] == "Healthy" for r in donors),
    analyzed_source_healthy_donors=sum(r["analyzed_donor"] and r["source_stage"] == "Healthy" for r in donors),
    analyzed_false_healthy=sum(r["analyzed_donor"] and r["condition_harmonized"] == "Healthy" and r["source_stage"] != "Healthy" for r in donors),
    analyzed_controls_labelled_masld=sum(r["analyzed_donor"] and r["condition_harmonized"] == "MASLD" and r["source_stage"] == "Healthy" for r in donors),
    table_s6_exact_join_donors=len(source_join),
    note="Donor identifiers are source-based, not genotype-verified. GSE189600 replicate independence is unresolved. No inference or expression refit.",
    environment=dict(python=sys.version, executable=sys.executable, platform=platform.platform(),
                     numpy=np.__version__, h5py=h5py.__version__, openpyxl=openpyxl.__version__),
)
with (OUT / "summary.json").open("x") as handle:
    json.dump(summary, handle, indent=2)
# Hash only small evidence inputs. Large atlas/count matrices are identified
# by path, size, mtime, and the header inspected, without reading their bodies.
inventory = []
for path in sorted(set(inputs + [atlas])):
    stat = path.stat()
    inventory.append(dict(path=str(path), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                          sha256=hashlib.sha256(path.read_bytes()).hexdigest() if stat.st_size < 5_000_000 else "not_hashed_large_input"))
write_table("input_inventory.tsv", inventory)
print(json.dumps(summary, indent=2))
