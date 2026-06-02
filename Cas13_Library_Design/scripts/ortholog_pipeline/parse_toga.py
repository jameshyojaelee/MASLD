#!/usr/bin/env python
"""
parse_toga.py — Parse TOGA Zoonomia mouse-human orthology output into the
standard L0_toga.tsv layer schema.

TOGA (Kirilenko et al. 2023, Science) provides whole-genome Cactus alignment +
ML-based orthology calls. The download for human(hg38)→mouse(mm39) consists of:
  - orthologsClassification.tsv.gz: t_gene, t_transcript, q_gene, q_transcript,
                                    orthology_class
  - geneAnnotation.bed.gz: bed12 with mouse mm39 coordinates per ortholog
                           projection (column 4 = ref_ENST.SYMBOL.chainID)
  - loss_summ_data.tsv.gz: per-entity intactness label (I/PI/UL/M/PM/L/PG)

TOGA's `q_gene` is an internal region ID (reg_NNN), NOT an Ensembl mouse ID.
We resolve the mouse Ensembl gene per projection by bedtools-intersecting the
mouse mm39 BED coordinates against mouse GENCODE vM38 gene annotation, taking
the gene with maximum overlap as the mouse ortholog locus.

Output schema (data/external/orthologs/layers/L0_toga.tsv):
  mouse_ensembl, mouse_symbol, human_ensembl, human_symbol,
  tier_H_toga, toga_orthology_class, confidence_tier, provenance_sources,
  notes

Confidence tier rules:
  one2one              → H  (gold-standard 1:1)
  one2many / many2one  → H  (paralog flag in notes)
  many2many            → H  (paralog flag in notes; multiple-row pair)
  one2zero             → not exported (no mouse ortholog)
  paralog_projection   → tier L w/ paralog flag (PG class only from loss_summ)
"""
from __future__ import annotations

import gzip
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
TOGA_DIR = PROJECT / "data/external/orthologs/curated_dbs/toga_zoonomia"
LAYERS = PROJECT / "data/external/orthologs/layers"
MOUSE_GTF = PROJECT / "data/ncrna_conservation/gencode.vM38.annotation.gtf.gz"
OUT_L0 = LAYERS / "L0_toga.tsv"

# Regex to parse q_transcript = <ref_ENST>.<ref_symbol>.<chain_id>
# The symbol may contain dots (e.g., MARCH1, MIR1-2HG, AC007326.13).
# We split on the first dot and the LAST dot; middle is the symbol.
TX_RE = re.compile(r"^([A-Z]+\d+(?:\.\d+)?)\.(.+)\.(\d+)$")


def parse_q_transcript(qt: str) -> tuple[str, str, str]:
    """Return (ref_enst, ref_symbol, chain_id). Symbol can contain dots."""
    m = TX_RE.match(qt)
    if not m:
        return ("", "", "")
    return m.group(1), m.group(2), m.group(3)


def build_mouse_gene_bed(gtf_path: Path, out_bed: Path) -> int:
    """Extract gene-level BED from mouse GENCODE GTF (chr, start, end, ensembl,
    symbol, strand, biotype, with stripped Ensembl version)."""
    n = 0
    with gzip.open(gtf_path, "rt") as f, open(out_bed, "w") as o:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, start, end, strand, attr = (
                fields[0], int(fields[3]) - 1, int(fields[4]), fields[6], fields[8]
            )
            gid_m = re.search(r'gene_id "([^"]+)"', attr)
            gn_m = re.search(r'gene_name "([^"]+)"', attr)
            gt_m = re.search(r'gene_type "([^"]+)"', attr)
            if not (gid_m and gn_m):
                continue
            gid = gid_m.group(1).split(".")[0]  # strip version
            sym = gn_m.group(1)
            biot = gt_m.group(1) if gt_m else ""
            o.write(f"{chrom}\t{start}\t{end}\t{gid}|{sym}|{biot}\t0\t{strand}\n")
            n += 1
    return n


def build_toga_projection_bed(toga_bed_gz: Path, out_bed: Path) -> int:
    """Extract a per-projection BED from TOGA's geneAnnotation.bed.gz.
    Keep cols 1..6, with name = q_transcript (matches orthologsClassification)."""
    n = 0
    with gzip.open(toga_bed_gz, "rt") as f, open(out_bed, "w") as o:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            o.write("\t".join(fields[:6]) + "\n")
            n += 1
    return n


def bedtools_intersect_best(toga_bed: Path, mouse_bed: Path,
                            tmp_dir: Path) -> pd.DataFrame:
    """Run bedtools intersect -wo to find mouse genes overlapping each TOGA
    projection, then keep the longest overlap per projection."""
    out = tmp_dir / "intersect.bed"
    # -s require same strand; without -s, TOGA may align on either strand
    cmd = ["bedtools", "intersect", "-a", str(toga_bed), "-b", str(mouse_bed),
           "-wo"]
    with open(out, "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)
    rows = []
    with open(out) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            # a_chrom,a_start,a_end,a_name,a_score,a_strand,
            # b_chrom,b_start,b_end,b_name,b_score,b_strand,overlap
            if len(f) < 13:
                continue
            q_tx = f[3]
            b_name = f[9]  # ensg|sym|biot
            overlap = int(f[12])
            parts = b_name.split("|", 2)
            if len(parts) != 3:
                continue
            ens, sym, biot = parts
            rows.append((q_tx, ens, sym, biot, overlap))
    if not rows:
        return pd.DataFrame(columns=["q_transcript", "mouse_ensembl",
                                     "mouse_symbol", "mouse_biotype",
                                     "overlap"])
    df = pd.DataFrame(rows, columns=["q_transcript", "mouse_ensembl",
                                     "mouse_symbol", "mouse_biotype",
                                     "overlap"])
    # Best (max overlap) ensembl per q_transcript
    df = df.sort_values(["q_transcript", "overlap"], ascending=[True, False])
    df = df.drop_duplicates(subset=["q_transcript"], keep="first")
    return df


def load_human_symbols() -> pd.DataFrame:
    """Build ENSG -> human_symbol map from existing biomaRt or gencode v49."""
    p = PROJECT / "data/gencode_v49_gene_metadata.tsv.gz"
    if p.exists():
        df = pd.read_csv(p, sep="\t", dtype=str, low_memory=False)
        rename = {}
        for c in df.columns:
            cl = c.lower()
            if cl in ("gene_id", "ensg", "human_ensembl", "ensembl_gene_id",
                       "gene_id_unversioned"):
                rename[c] = "human_ensembl"
            elif cl in ("gene_name", "symbol", "gene_symbol",
                         "human_symbol", "hgnc_symbol"):
                rename[c] = "human_symbol"
        df = df.rename(columns=rename)
        if "human_ensembl" not in df.columns:
            # try first column
            df = df.rename(columns={df.columns[0]: "human_ensembl"})
        df["human_ensembl"] = df["human_ensembl"].astype(str).str.split(".").str[0]
        if "human_symbol" not in df.columns:
            # second column as fallback symbol
            df["human_symbol"] = df.iloc[:, 1].astype(str)
        return df[["human_ensembl", "human_symbol"]].drop_duplicates(
            subset=["human_ensembl"], keep="first")
    # fallback to biomaRt L1
    p2 = PROJECT / "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz"
    df = pd.read_csv(p2, sep="\t", dtype=str)
    df = df.rename(columns={"human_gene_symbol": "human_symbol",
                            "human_gene_ensembl": "human_ensembl"})
    df["human_ensembl"] = df["human_ensembl"].astype(str).str.split(".").str[0]
    return df[["human_ensembl", "human_symbol"]].drop_duplicates(
        subset=["human_ensembl"], keep="first")


def main() -> int:
    if not TOGA_DIR.exists():
        print(f"[parse_toga] {TOGA_DIR} missing — run download first",
              file=sys.stderr)
        return 1

    ortho_path = TOGA_DIR / "orthologsClassification.tsv.gz"
    bed_path = TOGA_DIR / "geneAnnotation.bed.gz"

    print(f"[parse_toga] reading {ortho_path}...", file=sys.stderr)
    ortho = pd.read_csv(ortho_path, sep="\t", dtype=str)
    print(f"[parse_toga]   rows: {len(ortho):,}", file=sys.stderr)

    # Drop one2zero (no mouse ortholog)
    n0 = len(ortho)
    ortho = ortho[ortho["orthology_class"] != "one2zero"].copy()
    print(f"[parse_toga]   dropped {n0 - len(ortho):,} one2zero rows; "
          f"keeping {len(ortho):,}", file=sys.stderr)

    # Strip ENSG version
    ortho["human_ensembl"] = ortho["t_gene"].astype(str).str.split(".").str[0]

    # Extract embedded symbol + chain_id from q_transcript (used for human-side
    # symbol if reference annotation didn't provide one downstream)
    parts = ortho["q_transcript"].apply(parse_q_transcript)
    ortho["q_ref_enst"] = parts.apply(lambda t: t[0])
    ortho["q_ref_symbol"] = parts.apply(lambda t: t[1])  # human symbol embedded
    ortho["q_chain_id"] = parts.apply(lambda t: t[2])

    # Build mouse GENCODE BED + TOGA projection BED, intersect
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        mouse_bed = tmp_dir / "mouse_gencode_vM38_genes.bed"
        proj_bed = tmp_dir / "toga_projections.bed"
        print("[parse_toga] building mouse GENCODE vM38 gene BED...",
              file=sys.stderr)
        n_mouse = build_mouse_gene_bed(MOUSE_GTF, mouse_bed)
        print(f"[parse_toga]   {n_mouse:,} mouse genes in BED",
              file=sys.stderr)
        print("[parse_toga] building TOGA projection BED...", file=sys.stderr)
        n_proj = build_toga_projection_bed(bed_path, proj_bed)
        print(f"[parse_toga]   {n_proj:,} projection BED rows",
              file=sys.stderr)
        print("[parse_toga] bedtools intersect (best-overlap per projection)...",
              file=sys.stderr)
        proj_to_mouse = bedtools_intersect_best(proj_bed, mouse_bed, tmp_dir)
        print(f"[parse_toga]   resolved {len(proj_to_mouse):,} projections to "
              f"mouse Ensembl IDs", file=sys.stderr)

    # Merge mouse ENSMUSG onto orthology rows by q_transcript
    ortho = ortho.merge(proj_to_mouse, on="q_transcript", how="left")
    n_with_mouse = ortho["mouse_ensembl"].notna().sum()
    print(f"[parse_toga] projections with mouse ENSMUSG: "
          f"{n_with_mouse:,} / {len(ortho):,}", file=sys.stderr)

    # Per-(human_ensembl, mouse_ensembl) collapse: keep best (most-conservative)
    # class. Class priority: one2one > one2many ≈ many2one > many2many.
    class_rank = {"one2one": 0, "one2many": 1, "many2one": 1, "many2many": 2}
    ortho["_class_rank"] = ortho["orthology_class"].map(class_rank).fillna(99)
    # Drop rows lacking mouse_ensembl resolution
    pairs = ortho.dropna(subset=["mouse_ensembl"]).copy()
    pairs = pairs.sort_values(
        ["human_ensembl", "mouse_ensembl", "_class_rank"],
        ascending=[True, True, True]
    ).drop_duplicates(subset=["human_ensembl", "mouse_ensembl"], keep="first")

    # Add human symbol from GENCODE v49 metadata
    hsyms = load_human_symbols()
    pairs = pairs.merge(hsyms, on="human_ensembl", how="left")
    # Fall back to TOGA-embedded human symbol if GENCODE didn't have it
    pairs["human_symbol"] = pairs["human_symbol"].fillna(pairs["q_ref_symbol"])

    # Confidence tier + notes
    def _tier(cls):
        return "H"  # one2one/one2many/many2one/many2many all H given Cactus + ML

    def _notes(cls):
        if cls == "one2one":
            return ""
        if cls in ("one2many", "many2one"):
            return "paralog_branch"
        if cls == "many2many":
            return "paralog_many2many"
        return ""

    pairs["tier_H_toga"] = 1
    pairs["confidence_tier"] = pairs["orthology_class"].apply(_tier)
    pairs["provenance_sources"] = "toga_zoonomia"
    pairs["notes"] = pairs["orthology_class"].apply(_notes)
    pairs = pairs.rename(columns={"orthology_class": "toga_orthology_class"})

    out = pairs[[
        "mouse_ensembl", "mouse_symbol", "mouse_biotype",
        "human_ensembl", "human_symbol",
        "tier_H_toga", "toga_orthology_class", "confidence_tier",
        "provenance_sources", "notes",
    ]].copy()
    # Drop pairs missing either gene ID
    out = out.dropna(subset=["mouse_ensembl", "human_ensembl"])
    out = out.drop_duplicates(subset=["mouse_ensembl", "human_ensembl"],
                              keep="first")

    OUT_L0.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_L0, sep="\t", index=False)
    print(f"[parse_toga] wrote {OUT_L0}  ({len(out):,} pairs)",
          file=sys.stderr)
    print("[parse_toga] orthology class distribution:", file=sys.stderr)
    print(out["toga_orthology_class"].value_counts().to_string(),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
