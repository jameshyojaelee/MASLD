"""Build the symbol-based mouse↔human ortholog table.

Inputs
------
1. Archived Ensembl-only ortholog map (Ensembl Biomart export):
     archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz
     cols: mouse_ensembl_gene_id, human_ensembl_gene_id, orthology_type
2. Human GENCODE v49 gene metadata:
     data/gencode_v49_gene_metadata.tsv.gz
     cols: gene_id, gene_name, chromosome, gene_biotype, ensembl_base
3. Mouse GENCODE vM37 GTF (parsed for gene_id ↔ gene_name):
     ~/reference_genome/refdata-gex-GRCm39-2024-A/annotation/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf.gz
4. Cross-check / supplementary symbol pairs from the legacy mouse-human comparison:
     archive/lfc0.3_baseline_2026-04-27/integration_results/human_mouse_ortholog_comparison.csv

Output
------
data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz with columns:
  human_gene_symbol, human_gene_ensembl,
  mouse_gene_symbol, mouse_gene_ensembl,
  ortholog_type, conservation_score (NA — Ensembl Biomart export does not carry %ID)

Subagent: d5-ortholog-mapper.
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
HOME = Path("/gpfs/commons/home/jameslee")

ORTHO_TSV = PROJECT_ROOT / "archive/streamlit_deg_explorer/data/mouse_human_orthologs.tsv.gz"
HUMAN_META = PROJECT_ROOT / "data/gencode_v49_gene_metadata.tsv.gz"
MOUSE_GTF = HOME / "reference_genome/refdata-gex-GRCm39-2024-A/annotation/gencode.vM37.chr_patch_hapl_scaff.annotation.gtf.gz"
SUPP_PAIRS = PROJECT_ROOT / "archive/lfc0.3_baseline_2026-04-27/integration_results/human_mouse_ortholog_comparison.csv"

OUT_DIR = PROJECT_ROOT / "data/external/orthologs"
OUT_TSV = OUT_DIR / "mouse_human_orthologs_symbols.tsv.gz"


GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')
GENE_NAME_RE = re.compile(r'gene_name "([^"]+)"')
GENE_TYPE_RE = re.compile(r'gene_type "([^"]+)"')


def parse_mouse_gtf(path: Path) -> pd.DataFrame:
    """Stream the mouse GTF, keep one row per gene_id with its symbol + biotype."""
    seen: dict[str, tuple[str, str]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 8)
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            mid = GENE_ID_RE.search(attrs)
            mname = GENE_NAME_RE.search(attrs)
            mtype = GENE_TYPE_RE.search(attrs)
            if not mid or not mname:
                continue
            full = mid.group(1)
            base = full.split(".")[0]
            seen[base] = (mname.group(1), mtype.group(1) if mtype else "")
    df = pd.DataFrame(
        [(b, n, t) for b, (n, t) in seen.items()],
        columns=["mouse_ensembl", "mouse_symbol", "mouse_biotype"],
    )
    return df


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[1/5] Loading Ensembl ortholog map: {ORTHO_TSV}")
    ortho = pd.read_csv(ORTHO_TSV, sep="\t")
    print(f"      {len(ortho):,} ortholog rows")

    print(f"[2/5] Loading human GENCODE v49 metadata: {HUMAN_META}")
    hmeta = pd.read_csv(HUMAN_META, sep="\t")
    hmap = (
        hmeta[["ensembl_base", "gene_name", "gene_biotype"]]
        .drop_duplicates("ensembl_base")
        .rename(columns={
            "ensembl_base": "human_ensembl",
            "gene_name": "human_symbol",
            "gene_biotype": "human_biotype",
        })
    )
    print(f"      {len(hmap):,} unique human ensembl_base entries")

    print(f"[3/5] Parsing mouse GENCODE vM37 GTF for ENSMUSG → symbol")
    mmap = parse_mouse_gtf(MOUSE_GTF)
    print(f"      {len(mmap):,} unique mouse ensembl entries")

    print("[4/5] Joining tables")
    merged = (
        ortho.rename(columns={
            "mouse_ensembl_gene_id": "mouse_ensembl",
            "human_ensembl_gene_id": "human_ensembl",
            "orthology_type": "ortholog_type",
        })
        .merge(hmap, on="human_ensembl", how="left")
        .merge(mmap, on="mouse_ensembl", how="left")
    )

    # Try to recover missing mouse symbols from the supplementary comparison file
    if SUPP_PAIRS.exists():
        supp = pd.read_csv(SUPP_PAIRS, usecols=["mouse_gene_id", "mouse_symbol"]).dropna()
        supp = supp.rename(columns={"mouse_gene_id": "mouse_ensembl", "mouse_symbol": "mouse_symbol_supp"})
        supp = supp.drop_duplicates("mouse_ensembl")
        merged = merged.merge(supp, on="mouse_ensembl", how="left")
        before = merged["mouse_symbol"].isna().sum()
        merged["mouse_symbol"] = merged["mouse_symbol"].fillna(merged["mouse_symbol_supp"])
        after = merged["mouse_symbol"].isna().sum()
        print(f"      Mouse-symbol fillna via supp pairs: {before} → {after} NaN")
        merged = merged.drop(columns=["mouse_symbol_supp"])

    merged["conservation_score"] = pd.NA  # not available in this export
    out = merged[[
        "human_symbol", "human_ensembl",
        "mouse_symbol", "mouse_ensembl",
        "ortholog_type", "conservation_score",
    ]].rename(columns={
        "human_symbol": "human_gene_symbol",
        "human_ensembl": "human_gene_ensembl",
        "mouse_symbol": "mouse_gene_symbol",
        "mouse_ensembl": "mouse_gene_ensembl",
    })

    # Drop rows missing both Ensembl IDs (would be entirely unrecoverable)
    keep = out["human_gene_ensembl"].notna() & out["mouse_gene_ensembl"].notna()
    dropped = (~keep).sum()
    out = out[keep].copy()

    print(f"[5/5] Writing {OUT_TSV}")
    out.to_csv(OUT_TSV, sep="\t", index=False, compression="gzip")

    print("--- Summary ---")
    print(f"  rows written: {len(out):,}  (dropped {dropped} fully-NaN-Ensembl rows)")
    print(f"  ortholog_type breakdown:")
    print(out["ortholog_type"].value_counts(dropna=False).to_string())
    print(f"  rows with both symbols populated: "
          f"{(out['human_gene_symbol'].notna() & out['mouse_gene_symbol'].notna()).sum():,}")
    print(f"  rows missing human_symbol: {out['human_gene_symbol'].isna().sum():,}")
    print(f"  rows missing mouse_symbol: {out['mouse_gene_symbol'].isna().sum():,}")


if __name__ == "__main__":
    main()
