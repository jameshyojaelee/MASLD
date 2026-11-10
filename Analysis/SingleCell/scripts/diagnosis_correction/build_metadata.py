"""Build source-qualified RNA labels without modifying historical releases."""
import argparse
import collections
import csv
import hashlib
import json
import re
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[4]
REVIEW = ROOT / "Analysis/SingleCell/reviews/control-label-review-20261021"
OLD = ROOT / "Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014"


def read(path, sep="\t"):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=sep))


def write(path, rows, sep="\t"):
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter=sep)
        writer.writeheader()
        writer.writerows(rows)


def soft(dataset):
    result = {}
    text = (REVIEW / "sources" / f"{dataset}_family.soft").read_text()
    for block in text.split("^SAMPLE = ")[1:]:
        lines = block.splitlines()
        row = {"gsm": lines[0]}
        for line in lines:
            if line.startswith("!Sample_characteristics_ch1 = "):
                key, sep, value = line.split(" = ", 1)[1].partition(": ")
                if sep:
                    row[key.lower()] = value
        result[lines[0]] = row
    return result


def build(out):
    out.mkdir(parents=True, exist_ok=False)
    libs = read(REVIEW / "library_labels.tsv")
    donors = read(REVIEW / "donor_labels.tsv")
    prior = {r["donor"]: r for r in read(OLD / "metadata/current_source_diagnosis_donor_roster.tsv")}
    assert len(libs) == 275 and len(donors) == 104 and len(prior) == 102
    g202 = soft("GSE202379")
    patient = {}
    for row in g202.values():
        key = "GSE202379_" + row["patient id"]
        if key in patient:
            assert patient[key]["disease status"] == row["disease status"]
            assert patient[key]["saf score"] == row["saf score"]
        patient[key] = row
    wb = openpyxl.load_workbook(REVIEW / "sources/GSE244832_240605_18liver_metadata.xlsx", read_only=True, data_only=True)
    native244 = {r[4]: r for r in list(wb.active.values)[1:] if r[1]}
    wb.close()
    pair244 = read(ROOT / "data/GSE244832/metadata/donor_pairing.csv", ",")
    rna244 = {}
    corrected_pair = []
    for pair in pair244:
        native = native244[pair["rna_gsm"]]
        assert pair["rna_sample"] == native[1]
        rna244["GSE244832_" + pair["donor_id"]] = native
        corrected_pair.append(dict(pair, original_crosswalk_condition=pair["condition"],
                                   condition=native[3], rna_source_diagnosis=native[3],
                                   source_donor_id=native[2], pairing_identity_verified="FALSE"))
    etiology136 = soft("GSE136103")
    donor_libs = collections.defaultdict(list)
    for row in libs:
        donor_libs[row["donor"]].append(row)
    qualified = []
    for old in donors:
        donor, dataset = old["donor"], old["dataset"]
        stage = old["source_stage"]
        diagnosis = stage
        histology = ""
        etiology = "MASLD_source_cohort" if dataset in ("GSE202379", "GSE244832", "GSE174748", "GSE189600") else ""
        reason = "source_diagnosis"
        source = "GEO source study"
        source_id = ""
        if dataset == "GSE202379":
            native = patient[donor]
            diagnosis, histology = native["disease status"], native["saf score"]
            source_id = native["gsm"]
            source = "GSE202379 GEO patient diagnosis and SAF, retained separately"
            if donor == "GSE202379_P98":
                assert histology == "S1A1F1"
                stage, reason = "Unresolved", "Healthy_control_diagnosis_conflicts_with_SAF_S1A1F1"
            elif donor == "GSE202379_P30":
                assert histology == "S0A1F0"
                reason = "source_nonsteatotic_control_with_activity_A1"
        elif dataset == "GSE244832":
            native = rna244[donor]
            diagnosis, source_id = native[3], native[4]
            stage = {"NORMAL": "Healthy", "MASL": "Steatosis", "MASH": "Steatohepatitis"}[diagnosis]
            source = "GSE244832_240605_18liver_metadata.xlsx joined by RNA_ID and GSE_accession"
            reason = "deposited_RNA_diagnosis_independent_of_RNA_ATAC_crosswalk"
        elif dataset == "GSE136103":
            natives = [etiology136[r["sample"]] for r in donor_libs[donor]]
            causes = {n["cause of liver disease"] for n in natives}
            assert len(causes) == 1
            etiology = causes.pop()
            diagnosis = "Healthy" if stage == "Healthy" else etiology + " cirrhosis"
            source_id = ";".join(r["sample"] for r in donor_libs[donor])
            source = "GSE136103 GEO disease status and cause of liver disease"
        elif dataset == "GSE189600":
            stage = "Healthy" if old["condition"] == "Healthy" else "Steatohepatitis"
            diagnosis = "Healthy" if stage == "Healthy" else "NASH"
            reason = "replicate_independence_unresolved_sensitivity_only"
        elif dataset == "Liver_Atlas":
            diagnosis = stage
            histology = "steatosis=" + old["native_Steatosis"] + ";fibrosis=" + old["native_Fibrosis"]
            etiology = "unresolved"
            source = "Guilliams Cell 2022 Table S6 exact H-number join"
            source_id = donor.removeprefix("Liver_Atlas_")
            reason = "source_histology_only_all_19_excluded_from_disease_inference"
        if stage == "Healthy":
            binary = "Healthy"
        elif stage == "Unresolved":
            binary = "Unresolved"
        elif dataset == "Liver_Atlas" or (dataset == "GSE136103" and etiology != "NAFLD"):
            binary = "Other_liver_disease"
        else:
            binary = "MASLD"
        analyzed = old["analyzed_donor"] == "True"
        prior_exclude = prior.get(donor, {}).get("exclude_stage_analysis", "TRUE").lower() == "true"
        stage_exclude = prior_exclude or donor == "GSE202379_P98" or not analyzed
        binary_eligible = analyzed and binary in ("Healthy", "MASLD") and dataset not in ("Liver_Atlas", "GSE189600")
        row = dict(donor=donor, dataset=dataset, disease_stage_coarse=stage,
                   stage_ordinal={"Healthy": 0, "Steatosis": 1, "Steatohepatitis": 2, "Cirrhosis": 3}.get(stage, ""),
                   condition_harmonized=binary, source_diagnosis=diagnosis,
                   source_histology=histology, source_etiology=etiology,
                   disease_stage_source=source, source_accession=source_id,
                   qualification_reason=reason, analyzed_donor=analyzed,
                   exclude_stage_analysis=stage_exclude, exclude=stage_exclude,
                   binary_eligible=binary_eligible,
                   stage_primary_eligible=not stage_exclude and dataset != "GSE189600" and stage != "Cirrhosis",
                   strict_binary_eligible=binary_eligible and donor != "GSE202379_P30",
                   n_source_records=old["n_libraries"], n_cells=old["n_cells"])
        measured = re.fullmatch(r"S\d+A\d+F([0-4])", histology)
        row["F_stage_documented"] = int(measured.group(1)) if measured else ""
        qualified.append(row)
    by_donor = {r["donor"]: r for r in qualified}
    run_rows = []
    extended = {r["sample"]: r for r in read(OLD / "metadata/donor_metadata_stage_corrected.tsv")}
    extended_rows = []
    for lib in libs:
        q = by_donor[lib["donor"]]
        row = dict(lib)
        row.update(q)
        row["n_donor_cells"] = q["n_cells"]
        row["n_cells"] = lib["n_cells"]
        row["condition_original"] = lib["condition"]
        row["condition_harmonized_original"] = lib["condition_harmonized"]
        row["condition"] = {"Healthy": "Healthy", "Steatosis": "MASL", "Steatohepatitis": "MASH", "Cirrhosis": "Cirrhotic", "Unresolved": "Unresolved"}[q["disease_stage_coarse"]]
        run_rows.append(row)
        x = dict(extended.get(lib["sample"], {}))
        x.update(row)
        x["disease_stage_numeric"] = q["stage_ordinal"]
        x["condition_binary"] = q["condition_harmonized"]
        x["condition_binary_num"] = {"Healthy": 0, "MASLD": 1}.get(q["condition_harmonized"], "")
        if q["dataset"] == "GSE202379":
            x["F_stage_documented"] = q["F_stage_documented"]
            x["F_stage_doc_source"] = "deposited numeric SAF only; free-text healthy control is not F0"
        if lib["donor"] == "GSE202379_P98":
            x["F_stage_documented"] = 1
            x["F_stage_doc_source"] = "GSM6112244 measured SAF S1A1F1; diagnosis conflict retained"
        extended_rows.append(x)
    assert by_donor["GSE202379_P98"]["exclude_stage_analysis"]
    assert by_donor["GSE244832_D11"]["disease_stage_coarse"] == "Steatosis"
    assert by_donor["GSE244832_D14"]["disease_stage_coarse"] == "Steatohepatitis"
    assert all(not r["binary_eligible"] for r in qualified if r["source_etiology"] in ("Alcohol", "PBC"))
    write(out / "donor_metadata.tsv", qualified)
    write(out / "library_metadata.tsv", run_rows)
    write(out / "donor_metadata_stage_corrected.tsv", extended_rows)
    write(out / "current_source_diagnosis_donor_roster.tsv", [r for r in qualified if r["analyzed_donor"]])
    donor_compat = []
    for row in qualified:
        if row["analyzed_donor"]:
            x = dict(row)
            x["sample"] = x.pop("donor")
            donor_compat.append(x)
    write(out / "donor_metadata_stage_corrected_donor.tsv", donor_compat)
    write(out / "GSE244832_RNA_donor_pairing.csv", corrected_pair, ",")
    counts = collections.Counter((r["dataset"], r["disease_stage_coarse"]) for r in qualified if r["analyzed_donor"])
    write(out / "diagnosis_counts.tsv", [dict(dataset=k[0], stage=k[1], n_donors=v) for k, v in sorted(counts.items())])
    write(out / "source_files.tsv", [dict(path=str(p), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted((REVIEW / "sources").glob("*")) if p.is_file()])
    print(json.dumps(dict(libraries=len(run_rows), donors=len(qualified), analyzed=sum(r["analyzed_donor"] for r in qualified), stage_sensitivity=sum(not r["exclude_stage_analysis"] and r["disease_stage_coarse"] != "Cirrhosis" for r in qualified)), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    build(parser.parse_args().output)
