#!/usr/bin/env python
"""
parse_orthofinder.py — Parse OrthoFinder 3.x output into L8_orthofinder.tsv.

OrthoFinder 3.x (Emms & Kelly 2019 Genome Biol) output for 2-species runs:
  - Orthogroups/Orthogroups.txt — MCL-clustered orthogroups (all pairs)
  - Phylogenetic_Hierarchical_Orthogroups/N0.tsv — gene-tree-reconciled
    Hierarchical Orthogroups at the root node. These refine complex OGs
    (where a single OG contains paralogs) into true ortholog sub-groups.
  - Orthologues/*.tsv — may be empty for 2-species runs in v3.x

Strategy:
  1. Parse Orthogroups.txt for all cross-species OGs
  2. Parse N0.tsv for HOG-level refined ortholog groupings
  3. For OGs that appear in N0.tsv, use the HOG sub-groupings to produce
     refined ortholog pairs (gene-tree reconciliation separates paralogs)
  4. For OGs NOT in N0.tsv (simple 2-gene OGs), use direct pairing
  5. Classify orthology type (one2one/one2many/many2one/many2many)

Maps ENSG-based OrthoFinder IDs back to gene symbols and Ensembl IDs using
the ID maps produced by prep_orthofinder_fastas.py.

Output schema:
  mouse_ensembl, mouse_symbol, mouse_biotype, human_ensembl, human_symbol,
  human_biotype, tier_H_orthofinder, orthofinder_orthogroup,
  orthofinder_orthology_type, confidence_tier, provenance_sources, notes

Usage:
  python parse_orthofinder.py [--work-dir PATH]
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PIPE = (PROJECT / "Cas13_Library_Design/scripts/ortholog_pipeline")
DEFAULT_WORK = PIPE / "orthofinder_work"
LAYERS_OUT = PROJECT / "data/external/orthologs/layers"


def load_id_maps(work_dir: Path) -> tuple[dict, dict]:
    """Load human and mouse ID maps.

    Returns:
        human_map: {ensg_base: {symbol, ensp, enst, length}}
        mouse_map: {ensg_base: {symbol, ensp, enst, length}}
    """
    hmap = {}
    mmap = {}
    for path, target in [(work_dir / "id_maps/human_id_map.tsv", hmap),
                         (work_dir / "id_maps/mouse_id_map.tsv", mmap)]:
        if not path.exists():
            print(f"WARNING: ID map not found: {path}", file=sys.stderr)
            continue
        df = pd.read_csv(path, sep="\t", dtype=str)
        for _, row in df.iterrows():
            target[row["ensembl_gene"]] = {
                "symbol": row["gene_symbol"],
                "ensp": row["ensembl_protein"],
                "enst": row["ensembl_transcript"],
                "length": int(row["protein_length"]),
            }
    print(f"  Loaded {len(hmap)} human, {len(mmap)} mouse ID mappings",
          file=sys.stderr)
    return hmap, mmap


def find_results_dir(work_dir: Path) -> Path:
    """Find the OrthoFinder results directory (Results_MonDD or similar)."""
    results_base = work_dir / "results"
    if not results_base.exists():
        raise FileNotFoundError(f"No results directory at {results_base}")
    candidates = sorted(results_base.glob("Results_*"))
    if not candidates:
        if (results_base / "Orthogroups").exists():
            return results_base
        raise FileNotFoundError(
            f"No Results_* directory found in {results_base}. "
            f"Contents: {list(results_base.iterdir())}")
    return candidates[-1]


def extract_gene_id(of_id: str) -> str:
    """Extract the ENSG/ENSMUSG base from OrthoFinder's gene ID.

    Our headers are ENSG00000000003|TSPAN6 (no species prefix since we
    used -X flag implicitly via meaningful species names).
    """
    parts = of_id.strip()
    if "|" in parts:
        return parts.split("|")[0]
    return parts


def parse_orthogroups_txt(results_dir: Path
                          ) -> dict[str, dict[str, list[str]]]:
    """Parse Orthogroups.txt -> {OG_id: {human: [...], mouse: [...]}}.

    Orthogroups.txt format (OrthoFinder 3.x):
      OG0000000: ENSG...|SYM ENSMUSG...|SYM ...
    """
    og_path = results_dir / "Orthogroups" / "Orthogroups.txt"
    if not og_path.exists():
        raise FileNotFoundError(f"Orthogroups.txt not found at {og_path}")

    og_dict: dict[str, dict[str, list[str]]] = {}
    with open(og_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(": ", 1)
            og_id = parts[0]
            if len(parts) < 2:
                continue
            members = parts[1].split()
            human_genes = []
            mouse_genes = []
            for m in members:
                gid = extract_gene_id(m)
                if gid.startswith("ENSG"):
                    human_genes.append(gid)
                elif gid.startswith("ENSMUSG"):
                    mouse_genes.append(gid)
            og_dict[og_id] = {"human": human_genes, "mouse": mouse_genes}

    n_both = sum(1 for v in og_dict.values()
                 if v["human"] and v["mouse"])
    n_human = sum(1 for v in og_dict.values()
                  if v["human"] and not v["mouse"])
    n_mouse = sum(1 for v in og_dict.values()
                  if v["mouse"] and not v["human"])
    print(f"  Parsed {len(og_dict)} orthogroups: "
          f"{n_both} both species, {n_human} human-only, {n_mouse} mouse-only",
          file=sys.stderr)
    return og_dict


def parse_n0_hogs(results_dir: Path
                  ) -> dict[str, list[dict[str, list[str]]]]:
    """Parse N0.tsv -> {OG_id: [{human: [...], mouse: [...]}, ...]}.

    N0.tsv has gene-tree-reconciled HOGs. Multiple rows can share the
    same OG (different sub-HOGs within the same orthogroup). Each sub-HOG
    represents a refined ortholog group that gene-tree reconciliation
    determined to be a separate orthologous unit.

    Returns a dict mapping OG_id to a list of sub-HOG dicts.
    """
    n0_path = results_dir / "Phylogenetic_Hierarchical_Orthogroups" / "N0.tsv"
    if not n0_path.exists():
        print(f"  N0.tsv not found; using raw orthogroups only",
              file=sys.stderr)
        return {}

    df = pd.read_csv(n0_path, sep="\t", dtype=str)
    hog_dict: dict[str, list[dict[str, list[str]]]] = defaultdict(list)

    for _, row in df.iterrows():
        og_id = str(row["OG"]).strip()
        human_str = str(row.get("Human", "") or "")
        mouse_str = str(row.get("Mouse", "") or "")

        human_genes = []
        mouse_genes = []
        if human_str and human_str != "nan":
            human_genes = [extract_gene_id(g.strip())
                          for g in human_str.split(", ") if g.strip()]
        if mouse_str and mouse_str != "nan":
            mouse_genes = [extract_gene_id(g.strip())
                          for g in mouse_str.split(", ") if g.strip()]

        hog_dict[og_id].append({
            "hog": str(row.get("HOG", "")),
            "human": human_genes,
            "mouse": mouse_genes,
        })

    print(f"  Parsed {sum(len(v) for v in hog_dict.values())} HOG rows "
          f"across {len(hog_dict)} orthogroups", file=sys.stderr)
    return dict(hog_dict)


def classify_orthology_type(n_human: int, n_mouse: int) -> str:
    """Classify orthology based on count of human and mouse genes."""
    if n_human == 1 and n_mouse == 1:
        return "one2one"
    elif n_human == 1 and n_mouse > 1:
        return "one2many"
    elif n_human > 1 and n_mouse == 1:
        return "many2one"
    else:
        return "many2many"


def build_ortholog_table(results_dir: Path,
                        human_map: dict, mouse_map: dict) -> pd.DataFrame:
    """Build the L8 ortholog table from OrthoFinder output.

    Strategy:
    1. Parse Orthogroups.txt for all OGs
    2. Parse N0.tsv for HOG-level refinement
    3. For OGs in N0.tsv with cross-species HOGs, use HOG sub-groupings
       (gene-tree reconciliation resolves paralogs vs orthologs)
    4. For OGs NOT in N0.tsv, use direct pairing from Orthogroups.txt
    5. The orthology type is based on the sub-HOG counts (refined) or
       the OG counts (raw)
    """
    og_dict = parse_orthogroups_txt(results_dir)
    hog_dict = parse_n0_hogs(results_dir)

    records = []
    ogs_with_hog = 0
    ogs_without_hog = 0

    for og_id, members in og_dict.items():
        human_genes = members["human"]
        mouse_genes = members["mouse"]

        # Skip species-specific OGs
        if not human_genes or not mouse_genes:
            continue

        if og_id in hog_dict:
            # Use HOG sub-groupings for refined ortholog pairs
            ogs_with_hog += 1
            for sub_hog in hog_dict[og_id]:
                h_genes = sub_hog["human"]
                m_genes = sub_hog["mouse"]
                if not h_genes or not m_genes:
                    continue
                orth_type = classify_orthology_type(len(h_genes), len(m_genes))
                hog_id = sub_hog["hog"]
                for hid in h_genes:
                    hinfo = human_map.get(hid, {})
                    for mid in m_genes:
                        minfo = mouse_map.get(mid, {})
                        records.append({
                            "mouse_ensembl": mid,
                            "mouse_symbol": minfo.get("symbol", mid),
                            "mouse_biotype": "protein_coding",
                            "human_ensembl": hid,
                            "human_symbol": hinfo.get("symbol", hid),
                            "human_biotype": "protein_coding",
                            "tier_H_orthofinder": 1,
                            "orthofinder_orthogroup": og_id,
                            "orthofinder_hog": hog_id,
                            "orthofinder_orthology_type": orth_type,
                            "confidence_tier": "H",
                            "provenance_sources": "orthofinder_3.1.2",
                            "notes": "hog_refined",
                        })
        else:
            # Simple OG: no HOG refinement needed (trivial tree)
            ogs_without_hog += 1
            orth_type = classify_orthology_type(
                len(human_genes), len(mouse_genes))
            for hid in human_genes:
                hinfo = human_map.get(hid, {})
                for mid in mouse_genes:
                    minfo = mouse_map.get(mid, {})
                    records.append({
                        "mouse_ensembl": mid,
                        "mouse_symbol": minfo.get("symbol", mid),
                        "mouse_biotype": "protein_coding",
                        "human_ensembl": hid,
                        "human_symbol": hinfo.get("symbol", hid),
                        "human_biotype": "protein_coding",
                        "tier_H_orthofinder": 1,
                        "orthofinder_orthogroup": og_id,
                        "orthofinder_hog": "",
                        "orthofinder_orthology_type": orth_type,
                        "confidence_tier": "H",
                        "provenance_sources": "orthofinder_3.1.2",
                        "notes": "",
                    })

    print(f"  OGs with HOG refinement: {ogs_with_hog}", file=sys.stderr)
    print(f"  OGs without HOG (simple): {ogs_without_hog}", file=sys.stderr)

    df = pd.DataFrame(records)

    # Deduplicate
    before = len(df)
    df = df.drop_duplicates(
        subset=["mouse_ensembl", "human_ensembl"], keep="first"
    ).reset_index(drop=True)
    after = len(df)
    if before != after:
        print(f"  Dedup: {before} -> {after} pairs "
              f"({before - after} duplicates removed)", file=sys.stderr)

    return df


def compare_to_existing_layers(df: pd.DataFrame) -> dict:
    """Compare OrthoFinder pairs to existing biomaRt and TOGA layers.

    Returns a summary dict for the report.
    """
    print("\n=== Comparison with existing layers ===", file=sys.stderr)
    summary = {}

    of_pairs = set(zip(df["mouse_ensembl"], df["human_ensembl"]))
    summary["of_total"] = len(of_pairs)
    print(f"OrthoFinder L8: {len(of_pairs)} unique pairs", file=sys.stderr)

    of_121 = set(zip(
        df.loc[df["orthofinder_orthology_type"] == "one2one", "mouse_ensembl"],
        df.loc[df["orthofinder_orthology_type"] == "one2one", "human_ensembl"]))
    summary["of_121"] = len(of_121)

    # biomaRt
    bm_path = PROJECT / "data/external/orthologs/mouse_human_orthologs_symbols.tsv.gz"
    if bm_path.exists():
        bm = pd.read_csv(bm_path, sep="\t", dtype=str)
        bm["mouse_gene_ensembl"] = bm["mouse_gene_ensembl"].str.split(".").str[0]
        bm["human_gene_ensembl"] = bm["human_gene_ensembl"].str.split(".").str[0]
        bm_pairs = set(zip(bm["mouse_gene_ensembl"], bm["human_gene_ensembl"]))
        overlap = of_pairs & bm_pairs
        of_only = of_pairs - bm_pairs
        bm_only = bm_pairs - of_pairs
        bm_121 = set(zip(
            bm.loc[bm["ortholog_type"] == "ortholog_one2one", "mouse_gene_ensembl"],
            bm.loc[bm["ortholog_type"] == "ortholog_one2one", "human_gene_ensembl"]))
        overlap_121 = bm_121 & of_121
        jaccard = len(overlap) / len(of_pairs | bm_pairs) if of_pairs | bm_pairs else 0

        summary["bm_total"] = len(bm_pairs)
        summary["bm_overlap"] = len(overlap)
        summary["bm_of_only"] = len(of_only)
        summary["bm_only"] = len(bm_only)
        summary["bm_jaccard"] = jaccard
        summary["bm_121"] = len(bm_121)
        summary["bm_121_overlap"] = len(overlap_121)

        print(f"\nbiomaRt: {len(bm_pairs)} pairs", file=sys.stderr)
        print(f"  Overlap: {len(overlap)} ({len(overlap)/len(of_pairs)*100:.1f}% of OF)",
              file=sys.stderr)
        print(f"  OF-only: {len(of_only)}", file=sys.stderr)
        print(f"  biomaRt-only: {len(bm_only)}", file=sys.stderr)
        print(f"  Jaccard: {jaccard:.3f}", file=sys.stderr)
        print(f"  biomaRt 1:1: {len(bm_121)} | OF 1:1: {len(of_121)} | "
              f"1:1 overlap: {len(overlap_121)}", file=sys.stderr)

    # TOGA
    toga_path = LAYERS_OUT / "L0_toga.tsv"
    if toga_path.exists():
        tg = pd.read_csv(toga_path, sep="\t", dtype=str)
        tg["mouse_ensembl"] = tg["mouse_ensembl"].str.split(".").str[0]
        tg["human_ensembl"] = tg["human_ensembl"].str.split(".").str[0]
        tg_pairs = set(zip(tg["mouse_ensembl"], tg["human_ensembl"]))
        overlap_tg = of_pairs & tg_pairs
        of_only_tg = of_pairs - tg_pairs
        tg_only = tg_pairs - of_pairs
        jaccard_tg = len(overlap_tg) / len(of_pairs | tg_pairs) if of_pairs | tg_pairs else 0

        summary["toga_total"] = len(tg_pairs)
        summary["toga_overlap"] = len(overlap_tg)
        summary["toga_of_only"] = len(of_only_tg)
        summary["toga_only"] = len(tg_only)
        summary["toga_jaccard"] = jaccard_tg

        print(f"\nTOGA: {len(tg_pairs)} pairs", file=sys.stderr)
        print(f"  Overlap: {len(overlap_tg)} ({len(overlap_tg)/len(of_pairs)*100:.1f}% of OF)",
              file=sys.stderr)
        print(f"  OF-only: {len(of_only_tg)}", file=sys.stderr)
        print(f"  TOGA-only: {len(tg_only)}", file=sys.stderr)
        print(f"  Jaccard: {jaccard_tg:.3f}", file=sys.stderr)

    # BLAST RBH
    rbh_path = LAYERS_OUT / "L4_blast_rbh.tsv"
    if rbh_path.exists():
        rbh = pd.read_csv(rbh_path, sep="\t", dtype=str)
        rbh["mouse_ensembl"] = rbh["mouse_ensembl"].str.split(".").str[0]
        rbh["human_ensembl"] = rbh["human_ensembl"].str.split(".").str[0]
        rbh_pairs = set(zip(rbh["mouse_ensembl"], rbh["human_ensembl"]))
        overlap_rbh = of_pairs & rbh_pairs

        summary["rbh_total"] = len(rbh_pairs)
        summary["rbh_overlap"] = len(overlap_rbh)

        print(f"\nBLAST RBH: {len(rbh_pairs)} pairs", file=sys.stderr)
        print(f"  Overlap: {len(overlap_rbh)} ({len(overlap_rbh)/len(of_pairs)*100:.1f}% of OF)",
              file=sys.stderr)

    # Three-way overlap
    if bm_path.exists() and toga_path.exists():
        triple = of_pairs & bm_pairs & tg_pairs
        any_two = (of_pairs & bm_pairs) | (of_pairs & tg_pairs) | (bm_pairs & tg_pairs)
        summary["triple_overlap"] = len(triple)
        print(f"\nThree-way (OF + biomaRt + TOGA): {len(triple)}", file=sys.stderr)

    return summary


def spot_check_canonical_pairs(df: pd.DataFrame) -> list[dict]:
    """Verify canonical ortholog pairs are present and correctly classified."""
    print("\n=== Spot-check canonical pairs ===", file=sys.stderr)
    canonical = [
        ("TP53", "Trp53"),
        ("BRCA1", "Brca1"),
        ("ALB", "Alb"),
        ("PNPLA3", "Pnpla3"),
        ("HSD17B13", "Hsd17b13"),
        ("THRB", "Thrb"),
        ("HNF4A", "Hnf4a"),
        ("APOE", "Apoe"),
    ]
    results = []
    for hsym, msym in canonical:
        match = df[(df["human_symbol"] == hsym) & (df["mouse_symbol"] == msym)]
        if len(match) == 0:
            match = df[(df["human_symbol"].str.upper() == hsym.upper()) &
                      (df["mouse_symbol"].str.upper() == msym.upper())]
        if len(match) > 0:
            row = match.iloc[0]
            status = "OK"
            info = (f"OG={row['orthofinder_orthogroup']}, "
                    f"type={row['orthofinder_orthology_type']}")
            print(f"  {hsym}/{msym}: {info}  [{status}]", file=sys.stderr)
        else:
            status = "MISSING"
            info = "not found in OrthoFinder output"
            print(f"  {hsym}/{msym}: NOT FOUND  [WARN]", file=sys.stderr)
        results.append({"human": hsym, "mouse": msym,
                       "status": status, "info": info})
    return results


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work-dir", type=Path, default=DEFAULT_WORK,
                   help="OrthoFinder work directory")
    p.add_argument("--results-dir", type=Path, default=None,
                   help="Override: direct path to OrthoFinder Results_* dir")
    p.add_argument("--out", type=Path,
                   default=LAYERS_OUT / "L8_orthofinder.tsv",
                   help="Output TSV path")
    args = p.parse_args()

    print("=== parse_orthofinder.py ===", file=sys.stderr)

    # Load ID maps
    print("\n1. Loading ID maps...", file=sys.stderr)
    human_map, mouse_map = load_id_maps(args.work_dir)

    # Find results directory
    print("\n2. Finding OrthoFinder results...", file=sys.stderr)
    if args.results_dir:
        rdir = args.results_dir
    else:
        rdir = find_results_dir(args.work_dir)
    print(f"  Results dir: {rdir}", file=sys.stderr)

    # Build ortholog table
    print("\n3. Building ortholog table...", file=sys.stderr)
    df = build_ortholog_table(rdir, human_map, mouse_map)

    # Summary statistics
    print(f"\n=== Summary ===", file=sys.stderr)
    print(f"Total pairs: {len(df):,}", file=sys.stderr)
    print(f"Unique human genes: {df['human_ensembl'].nunique():,}",
          file=sys.stderr)
    print(f"Unique mouse genes: {df['mouse_ensembl'].nunique():,}",
          file=sys.stderr)
    print(f"\nOrthology type distribution:", file=sys.stderr)
    print(df["orthofinder_orthology_type"].value_counts().to_string(),
          file=sys.stderr)
    print(f"\nOrthogroups: {df['orthofinder_orthogroup'].nunique():,}",
          file=sys.stderr)

    # Spot-check canonical pairs
    spot_results = spot_check_canonical_pairs(df)

    # Compare to existing layers
    comp = compare_to_existing_layers(df)

    # Write output
    # Drop the orthofinder_hog column (internal-only)
    out_cols = ["mouse_ensembl", "mouse_symbol", "mouse_biotype",
                "human_ensembl", "human_symbol", "human_biotype",
                "tier_H_orthofinder", "orthofinder_orthogroup",
                "orthofinder_orthology_type", "confidence_tier",
                "provenance_sources", "notes"]
    out_df = df[[c for c in out_cols if c in df.columns]]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"\nWrote {len(out_df):,} pairs to {args.out}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
