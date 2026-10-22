"""Compare existing diagnoses to freshly downloaded primary-source records."""
import csv
import hashlib
import json
from pathlib import Path

import openpyxl

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[3]


def table(path, sep="\t"):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=sep))


def write(name, rows):
    with (OUT / name).open("x", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def soft_samples(dataset):
    text = (OUT / "sources" / f"{dataset}_family.soft").read_text()
    result = {}
    for block in text.split("^SAMPLE = ")[1:]:
        lines = block.splitlines()
        record = {"gsm": lines[0]}
        for line in lines:
            if line.startswith("!Sample_title = "):
                record["title"] = line.split(" = ", 1)[1]
            elif line.startswith("!Sample_characteristics_ch1 = "):
                key, sep, value = line.split(" = ", 1)[1].partition(": ")
                if sep:
                    record[key.lower()] = value
        result[record["gsm"]] = record
    return result


donors = {r["donor"]: r for r in table(OUT / "donor_labels.tsv")}
g202 = soft_samples("GSE202379")
g244 = soft_samples("GSE244832")
run_gsm = {}
# The cached SraRunInfo.csv is headerless; parse positions, retaining row 1.
with (ROOT / "data/GSE202379/metadata/SraRunInfo.csv").open() as handle:
    for row in csv.reader(handle):
        assert row[0].startswith("SRR") and row[11].startswith("GSM")
        run_gsm[row[0]] = row[11]

rows = []
for pairing in table(ROOT / "data/GSE202379/metadata/donor_pairing.csv", ","):
    donor = "GSE202379_" + pairing["donor_id"]
    for run in pairing["rna_srrs"].split(";"):
        record = g202[run_gsm[run]]
        assert record["patient id"] == pairing["donor_id"]
        assert record["disease status"] == pairing["condition"]
        saf = record["saf score"]
        rows.append(dict(dataset="GSE202379", donor=donor, library=run,
                         gsm=record["gsm"], source_title=record["title"],
                         deposited_diagnosis=record["disease status"],
                         deposited_histology=saf, current_stage=donors[donor]["source_stage"],
                         finding=("control_label_conflicts_with_SAF_S1A1F1" if pairing["donor_id"] == "P98" else
                                  "nonsteatotic_control_has_activity_A1" if pairing["donor_id"] == "P30" else "diagnosis_join_agrees")))

path = OUT / "sources/GSE244832_240605_18liver_metadata.xlsx"
workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
native = {}
for row in list(workbook.active.values)[1:]:
    if row[1]:
        assert row[1] not in native
        native[row[1]] = row
workbook.close()
assert len(native) == 18
diagnoses = []
for pairing in table(ROOT / "data/GSE244832/metadata/donor_pairing.csv", ","):
    record = native[pairing["rna_sample"]]
    assert record[4] == pairing["rna_gsm"]
    assert g244[record[4]]["title"].split(",")[0] == pairing["rna_sample"]
    donor = "GSE244832_" + pairing["donor_id"]
    agrees = pairing["condition"] == record[3]
    diagnoses.append(dict(donor=donor, rna_id=record[1], gsm=record[4],
                           public_donor_id=record[2], local_diagnosis=pairing["condition"],
                           deposited_diagnosis=record[3], agrees=agrees))
    for run in pairing["rna_srrs"].split(";"):
        rows.append(dict(dataset="GSE244832", donor=donor, library=run,
                         gsm=record[4], source_title=g244[record[4]]["title"],
                         deposited_diagnosis=record[3], deposited_histology="not_in_GEO_spreadsheet",
                         current_stage=donors[donor]["source_stage"],
                         finding="diagnosis_join_agrees" if agrees else "local_diagnosis_disagrees_with_GEO_supplement"))
write("GSE244832_native_diagnosis_comparison.tsv", diagnoses)
write("primary_source_diagnosis_checks.tsv", rows)

# Existing donor-deletion results show actual fitted membership in H59; these
# rows are inspection of stored results, not a newly run sensitivity analysis.
deletions = table(ROOT / "Analysis/Hypotheses/H59-20261015T112432Z/out/outcome/E7_single_donor_deletions.tsv")
selected = [r for r in deletions if r["deleted_donor"] in (
    "GSE202379_P98", "GSE202379_P30", "GSE244832_D11", "GSE244832_D14")]
write("existing_H59_deletion_evidence.tsv", selected)
source_files = sorted((OUT / "sources").glob("*"))
write("downloaded_source_checksums.tsv", [dict(file=str(p.relative_to(OUT)), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in source_files])
print(json.dumps(dict(gse202379_checked_libraries=sum(r["dataset"] == "GSE202379" for r in rows),
                      gse244832_checked_libraries=sum(r["dataset"] == "GSE244832" for r in rows),
                      gse244832_diagnosis_conflicts=[r for r in diagnoses if not r["agrees"]],
                      existing_H59_rows_for_flagged_donors=len(selected)), indent=2))
