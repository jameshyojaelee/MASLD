#!/usr/bin/env python3
"""Inventory Yakubovsky MATLAB objects without loading their large arrays."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import pandas as pd


TOKENS = {
    "expression": ("expr", "count", "gene", "rna"),
    "barcode": ("barcode", "spot"),
    "sample": ("sample", "donor", "patient"),
    "zonation": ("zon", "portal", "central", "periportal", "pericentral"),
    "lipid": ("lipid", "fat", "droplet", "steatos"),
}


def inventory_hdf5(path: Path) -> tuple[pd.DataFrame, dict[str, bool]]:
    rows: list[dict[str, object]] = []
    hits = {key: False for key in TOKENS}
    with h5py.File(path, "r") as handle:
        def visitor(name: str, obj: h5py.Dataset | h5py.Group) -> None:
            lower = name.lower()
            for key, words in TOKENS.items():
                hits[key] = hits[key] or any(word in lower for word in words)
            if isinstance(obj, h5py.Dataset):
                rows.append(
                    {
                        "file": path.name,
                        "path": name,
                        "object_type": "dataset",
                        "shape": "x".join(str(value) for value in obj.shape),
                        "dtype": str(obj.dtype),
                        "storage_bytes": int(obj.id.get_storage_size()),
                    }
                )
            else:
                rows.append(
                    {
                        "file": path.name,
                        "path": name,
                        "object_type": "group",
                        "shape": "",
                        "dtype": "",
                        "storage_bytes": 0,
                    }
                )
        handle.visititems(visitor)
    return pd.DataFrame(rows), hits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    inventories: list[pd.DataFrame] = []
    token_rows: list[dict[str, object]] = []
    for filename in ("v.mat", "zon_struct_all_full.mat"):
        path = args.data_dir / filename
        try:
            inventory, hits = inventory_hdf5(path)
            inventories.append(inventory)
            fmt = "matlab_v7.3_hdf5"
        except OSError as error:
            hits = {key: False for key in TOKENS}
            fmt = f"not_hdf5:{type(error).__name__}"
        for category, present in hits.items():
            token_rows.append(
                {
                    "file": filename,
                    "format": fmt,
                    "category": category,
                    "path_name_hit": present,
                }
            )

    if inventories:
        pd.concat(inventories, ignore_index=True).to_csv(
            args.output_dir / "hdf5_inventory.tsv", sep="\t", index=False
        )
    pd.DataFrame(token_rows).to_csv(
        args.output_dir / "schema_token_audit.tsv", sep="\t", index=False
    )

    metadata_path = args.data_dir / "human_samples_metadata.xlsx"
    metadata = pd.read_excel(metadata_path, sheet_name=None)
    metadata_inventory = []
    for sheet, frame in metadata.items():
        metadata_inventory.append(
            {
                "sheet": sheet,
                "n_rows": len(frame),
                "n_columns": frame.shape[1],
                "columns": "|".join(str(column) for column in frame.columns),
            }
        )
        frame.to_csv(
            args.output_dir / f"metadata_{sheet.replace(' ', '_')}.tsv",
            sep="\t",
            index=False,
        )
    pd.DataFrame(metadata_inventory).to_csv(
        args.output_dir / "metadata_inventory.tsv", sep="\t", index=False
    )

    all_hits = pd.DataFrame(token_rows).groupby("category")["path_name_hit"].any()
    preflight = {
        "status": "schema_inventory_complete_manual_join_mapping_required",
        "expression_name_hit": bool(all_hits.get("expression", False)),
        "barcode_name_hit": bool(all_hits.get("barcode", False)),
        "sample_name_hit": bool(all_hits.get("sample", False)),
        "zonation_name_hit": bool(all_hits.get("zonation", False)),
        "lipid_name_hit": bool(all_hits.get("lipid", False)),
        "authoritative_barcode_join_verified": False,
        "source_gate_pass": False,
        "next_action": "map MATLAB references and reconcile paper lipid donors before reading outcome arrays",
    }
    (args.output_dir / "schema_preflight.json").write_text(
        json.dumps(preflight, indent=2) + "\n"
    )
    pd.DataFrame([preflight]).to_csv(
        args.output_dir / "gate_status.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
