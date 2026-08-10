#!/usr/bin/env python3
"""Audit GSE287826 for an explicit AOI-to-donor mapping; never infer pairs."""

from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd


DONOR_RE = re.compile(r"(?:donor|patient|subject|individual)(?:[_ .-]*id)?", re.I)
NTC_RE = re.compile(r"(?:no[- ]?template|\bntc\b)", re.I)
HEALTHY_RE = re.compile(r"\bhealthy\b", re.I)
MASH_RE = re.compile(r"\b(?:mash|nash)\b", re.I)


def parse_soft(path: Path) -> pd.DataFrame:
    samples: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    characteristics: list[str] = []
    descriptions: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["characteristics"] = " | ".join(characteristics)
                    current["description"] = " | ".join(descriptions)
                    samples.append(current)
                current = {"sample_record": line.split("=", 1)[1].strip()}
                characteristics = []
                descriptions = []
            elif current is not None and line.startswith("!Sample_geo_accession = "):
                current["geo_accession"] = line.split("=", 1)[1].strip()
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                characteristics.append(line.split("=", 1)[1].strip())
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split("=", 1)[1].strip())
    if current is not None:
        current["characteristics"] = " | ".join(characteristics)
        current["description"] = " | ".join(descriptions)
        samples.append(current)

    frame = pd.DataFrame(samples).fillna("")
    if frame.empty:
        return frame
    combined = (
        frame.get("title", "")
        + " | "
        + frame.get("characteristics", "")
        + " | "
        + frame.get("description", "")
    )
    frame["is_ntc"] = combined.str.contains(NTC_RE)
    frame["condition"] = "unknown"
    frame.loc[combined.str.contains(HEALTHY_RE), "condition"] = "Healthy"
    frame.loc[combined.str.contains(MASH_RE), "condition"] = "MASH"
    return frame


def workbook_donor_candidates(sheets: dict[str, pd.DataFrame]) -> list[dict[str, object]]:
    candidates: list[dict[str, object]] = []
    for sheet, frame in sheets.items():
        for column in frame.columns:
            if DONOR_RE.search(str(column)):
                values = frame[column].dropna().astype(str)
                candidates.append(
                    {
                        "source": f"workbook:{sheet}",
                        "field": str(column),
                        "orientation": "column",
                        "n_values": len(values),
                        "n_unique": values.nunique(),
                        "values": "|".join(values.tolist()),
                    }
                )
        if frame.shape[1] > 1:
            first = frame.iloc[:, 0].astype(str)
            for row_index in frame.index[first.str.contains(DONOR_RE, na=False)]:
                values = frame.loc[row_index].iloc[1:].dropna().astype(str)
                candidates.append(
                    {
                        "source": f"workbook:{sheet}",
                        "field": str(frame.loc[row_index].iloc[0]),
                        "orientation": "row",
                        "n_values": len(values),
                        "n_unique": values.nunique(),
                        "values": "|".join(values.tolist()),
                    }
                )
    return candidates


def soft_donor_candidates(soft: pd.DataFrame) -> list[dict[str, object]]:
    if soft.empty:
        return []
    key_values: dict[str, list[str]] = {}
    for characteristics in soft["characteristics"]:
        for item in str(characteristics).split(" | "):
            if ":" not in item:
                continue
            key, value = item.split(":", 1)
            if DONOR_RE.search(key):
                key_values.setdefault(key.strip(), []).append(value.strip())
    return [
        {
            "source": "GEO_SOFT",
            "field": key,
            "orientation": "sample_characteristic",
            "n_values": len(values),
            "n_unique": len(set(values)),
            "values": "|".join(values),
        }
        for key, values in key_values.items()
    ]


def exact_gate(candidates: list[dict[str, object]]) -> tuple[bool, str]:
    for candidate in candidates:
        values = [value for value in str(candidate["values"]).split("|") if value]
        counts = Counter(values)
        if len(values) == 30 and len(counts) == 15 and set(counts.values()) == {2}:
            return True, f"explicit_30_to_15_mapping:{candidate['source']}:{candidate['field']}"
    return False, "no_explicit_30_to_15_two_aoi_donor_field"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", required=True, type=Path)
    parser.add_argument("--soft", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sheets = pd.read_excel(args.workbook, sheet_name=None)
    inventory = []
    for name, frame in sheets.items():
        inventory.append(
            {
                "sheet": name,
                "n_rows": len(frame),
                "n_columns": frame.shape[1],
                "columns": "|".join(str(column) for column in frame.columns),
            }
        )
    pd.DataFrame(inventory).to_csv(
        args.output_dir / "workbook_inventory.tsv", sep="\t", index=False
    )

    soft = parse_soft(args.soft)
    soft.to_csv(args.output_dir / "soft_sample_metadata.tsv", sep="\t", index=False)
    candidates = workbook_donor_candidates(sheets) + soft_donor_candidates(soft)
    candidate_columns = [
        "source", "field", "orientation", "n_values", "n_unique", "values"
    ]
    pd.DataFrame(candidates, columns=candidate_columns).to_csv(
        args.output_dir / "donor_key_candidates.tsv", sep="\t", index=False
    )

    n_soft = len(soft)
    n_ntc = int(soft["is_ntc"].sum()) if not soft.empty else 0
    biological = soft.loc[~soft["is_ntc"]] if not soft.empty else soft
    condition_counts = biological["condition"].value_counts().to_dict() if not soft.empty else {}
    donor_pass, reason = exact_gate(candidates)
    structure_pass = (
        n_soft == 31
        and n_ntc == 1
        and len(biological) == 30
        and condition_counts.get("Healthy", 0) == 8
        and condition_counts.get("MASH", 0) == 22
    )
    final_pass = structure_pass and donor_pass
    status = "pass" if final_pass else "skipped_no_donor_key" if structure_pass else "fail_source_structure"

    gate = {
        "dataset_id": "GSE287826",
        "status": status,
        "structure_pass": structure_pass,
        "authoritative_donor_key_pass": donor_pass,
        "reason": reason,
        "n_geo_records": n_soft,
        "n_ntc": n_ntc,
        "n_biological_aois": len(biological),
        "n_healthy_aois": condition_counts.get("Healthy", 0),
        "n_mash_aois": condition_counts.get("MASH", 0),
        "inference_authorized": final_pass,
        "pairing_inference_prohibited": True,
    }
    pd.DataFrame([gate]).to_csv(args.output_dir / "gate_status.tsv", sep="\t", index=False)
    (args.output_dir / "gate_status.json").write_text(json.dumps(gate, indent=2) + "\n")


if __name__ == "__main__":
    main()
