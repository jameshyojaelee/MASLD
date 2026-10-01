"""Save person/file joins from selected metadata fields and protein headers only."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / "data/PXD051911"
META_SHA = "340874d4f720f909fcfca428fed6da16ba3f7cfdf4d9fd53dd6e0378198ba2bc"
FILES = ("liver_protein_quant.txt", "plasma_protein_quant.txt", "plasma_protein_quant_filtered_and_batch_corrected.txt")
FIELDS = ["unique_identifier", "patient_name", "sample_group", "liver_proteomics_filename",
          "plasma_proteomics_filename", "plasma_batch_effects"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, records, fields):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if sha(SOURCE/"meta_data.txt") != META_SHA:
        raise ValueError("Fixed source identity metadata changed")
    headers, receipts, column_records = {}, {}, []
    for name in FILES:
        path = SOURCE/name
        with path.open("rb") as handle:
            raw = handle.readline(100001)
        if len(raw) > 100000 or not raw.endswith(b"\n"):
            raise ValueError("Header byte cap or termination invalid")
        columns = raw.decode().rstrip("\r\n").split("\t")
        expected = ["ProteinAccessions"] if "corrected" in name else ["ProteinAccessions", "Genes", "ProteinDescriptions"]
        if columns[:len(expected)] != expected or len(columns) != len(set(columns)):
            raise ValueError("Source feature/sample header differs or duplicated")
        headers[name] = set(columns[len(expected):])
        receipts[name] = dict(path=str(path), bytes_on_disk=path.stat().st_size,
                             header_bytes=len(raw), header_sha256=hashlib.sha256(raw).hexdigest(),
                             feature_annotation_columns=expected, sample_columns=len(headers[name]),
                             molecular_rows_read=False)
        column_records.extend(dict(source_file=name, column_index=i, column_name=c,
                                   role="feature_annotation" if i<len(expected) else "sample_identifier") for i,c in enumerate(columns))
    # Installed pandas usecols prevents clinical fields from being materialized.
    import pandas as pd
    meta = pd.read_csv(SOURCE/"meta_data.txt", sep="\t", usecols=FIELDS, dtype=str, keep_default_na=False)
    paired = []
    for r in meta.to_dict("records"):
        if r["sample_group"] != "initial_sample":
            continue
        if r["liver_proteomics_filename"] not in headers[FILES[0]] or r["plasma_proteomics_filename"] not in headers[FILES[1]]:
            continue
        paired.append(dict(**r, liver_column=r["liver_proteomics_filename"], plasma_column=r["plasma_proteomics_filename"],
                           corrected_plasma_column=r["unique_identifier"] if r["unique_identifier"] in headers[FILES[2]] else "",
                           liver_batch_title="2019" if r["liver_proteomics_filename"].startswith("2019") else
                           "2020" if r["liver_proteomics_filename"].startswith("2020") else "unresolved"))
    if len(paired)!=41:
        raise ValueError("Source header-verified initial pair count differs from41")
    for key in ("unique_identifier", "patient_name", "liver_column", "plasma_column"):
        if len({r[key] for r in paired}) != 41:
            raise ValueError("Ambiguous pair ownership: "+key)
    args.output.mkdir(parents=True, exist_ok=False)
    write(args.output/"codex_plasma_liver_source_contract_pairs.tsv", sorted(paired,key=lambda r:r["patient_name"]), list(paired[0]))
    write(args.output/"codex_plasma_liver_source_contract_headers.tsv", column_records, list(column_records[0]))
    producers = [ROOT/"scripts/analysis/histology_anchored_continuum/translation/02_paired_protein.R",
                 ROOT/"Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R",
                 ROOT/"Analysis/MASLD_Model_Benchmark/scripts/build_pxd051911_protein_transport_fixture.py",
                 ROOT/"scripts/manuscript/program_observability_map/t2_proteome_axis_map.R"]
    receipt = dict(paired_initial_people=41, same_person_same_source_visit=True, source_metadata_sha256=META_SHA,
                   headers=receipts, script_sha256=sha(Path(__file__)), producer_sha256={str(p):sha(p) for p in producers},
                   molecular_rows_read=False, clinical_values_selected=False, feature_row_axis_extracted=False,
                   intensity_missingness_assessed=False, independent_validation_asserted=False,
                   outputs_sha256={p.name:sha(p) for p in args.output.iterdir()})
    (args.output/"codex_plasma_liver_source_contract_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(dict(paired_initial_people=41, corrected_plasma_pair_columns=sum(bool(r["corrected_plasma_column"]) for r in paired),
                         sample_columns={n:r["sample_columns"] for n,r in receipts.items()}),sort_keys=True))


if __name__ == "__main__":
    main()
