#!/usr/bin/env python
"""02a: export the feature names each non-R assay actually measures.

The universe is an intersection, so it can only be built once every assay has
said what it measures. Spatial objects live in h5ad and are read once here, in
the spatial environment, and written as plain text so the R steps never need
that environment. Reading `var_names` in backed mode keeps a 10 GB object off
the node's memory.
"""
import argparse
import json
import pathlib
import sys

import anndata as ad

ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
CONTRACT = json.loads((ROOT / "scripts/analysis/cross_assay_modules/00_contract.json").read_text())


def gene_names(h5ad_path: pathlib.Path) -> list:
    a = ad.read_h5ad(h5ad_path, backed="r")
    names = [str(v) for v in a.var_names]
    try:
        a.file.close()
    except Exception:
        pass
    return names


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", required=True, type=pathlib.Path)
    args = ap.parse_args()
    out = args.out_root / "universe"
    out.mkdir(parents=True, exist_ok=True)

    summary = {}
    for key, label in (("visium_guilliams", "visium_guilliams"),
                       ("visium_vu", "visium_vu"),
                       ("cosmx", "cosmx_govaere")):
        path = ROOT / CONTRACT["inputs"][key]
        names = gene_names(path)
        target = out / f"measured_{label}.txt"
        if target.exists():
            raise SystemExit(f"Refusing to overwrite {target}")
        target.write_text("\n".join(names) + "\n")
        summary[label] = {"n_features": len(names), "source": str(path)}
        print(f"[02a] {label}: {len(names)} features", flush=True)

    (out / "measured_spatial_summary.json").write_text(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
