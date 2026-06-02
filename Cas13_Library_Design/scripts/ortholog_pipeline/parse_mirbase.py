#!/usr/bin/env python
"""
parse_mirbase.py — Build L2 miRNA ortholog layer from miRBase 22.1.

miRBase naming convention:
  hsa-miR-122-5p ↔ mmu-miR-122-5p   (same family, same arm)
  hsa-let-7a-5p  ↔ mmu-let-7a-5p
  hsa-let-7a-2-3p ↔ mmu-let-7a-2-3p (paralog cluster 2 of let-7a, 3p arm)

Strategy:
  1. Read miRBase mature.fa
  2. Extract miRNA name (e.g., "hsa-miR-122-5p")
  3. Split into stem family + optional paralog cluster id + arm:
     "hsa-let-7a-2-3p" -> family "let-7a", paralog_id "2", arm "3p"
     "hsa-miR-122-5p"  -> family "miR-122", paralog_id "", arm "5p"
     "hsa-miR-30c-2-3p"-> family "miR-30c", paralog_id "2", arm "3p"
  4. For each (family, paralog_id, arm) combo, pair hsa <-> mmu
  5. Map mature names to GENCODE miRNA gene names. Paralog-bearing rows ALSO
     check single-form symbols (paralogs share the same mature sequence at the
     family-arm level; HGNC concatenates the paralog number, MGI hyphenates).
  6. Collapse 5p / 3p / unknown arms at the (mouse_ensembl, human_ensembl)
     gene-pair level and emit common-schema TSV.

This is Tier H (miRBase family ID matching is the gold standard for miRNA
orthology — same family ID = orthologous mature miRNA across species).
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

PROJECT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
DEFAULT_FA = PROJECT / "data/external/orthologs/curated_dbs/mirbase_mature.fa.gz"
DEFAULT_OUT = PROJECT / "data/external/orthologs/layers/L2_mirbase.tsv"
HUMAN_GTF_META = PROJECT / "data/gencode_v49_gene_metadata.tsv.gz"
MOUSE_GTF_META = PROJECT / "Cas13_Library_Design/data/mouse_gencode_vM38_gene_metadata.csv"


def split_family_paralog_arm(rest: str) -> tuple[str, str, str]:
    """Split miRBase suffix into (family, paralog_id, arm).

    Inputs (the part after the species prefix has already been stripped):
      "miR-122-5p"      -> ("miR-122", "",  "5p")
      "let-7a-5p"       -> ("let-7a",  "",  "5p")
      "let-7a-2-3p"     -> ("let-7a",  "2", "3p")
      "miR-30c-2-3p"    -> ("miR-30c", "2", "3p")
      "miR-19b-1-5p"    -> ("miR-19b", "1", "5p")
      "miR-16-1-3p"     -> ("miR-16",  "1", "3p")    (numeric family + paralog)
      "miR-7-1-5p"      -> ("miR-7",   "1", "5p")
      "miR-1306"        -> ("miR-1306","",  "unknown")    (no arm, no paralog)
      "let-7c"          -> ("let-7c",  "",  "unknown")
    Rules:
      - Trailing "-5p" or "-3p" -> arm
      - After stripping arm, split on "-". A trailing pure-digit segment is the
        paralog cluster id ONLY when there are >=3 segments total
        (genus + family-id + paralog). For 2-segment stems like "miR-122",
        the digits are the family id itself, NOT a paralog.
    """
    arm_match = re.search(r"-(5p|3p)$", rest)
    if arm_match:
        arm = arm_match.group(1)
        stem = rest[: -(len(arm) + 1)]  # drop "-5p" / "-3p"
    else:
        arm = "unknown"
        stem = rest
    segs = stem.split("-")
    # Paralog id only if we have >=3 segments AND last segment is pure digits.
    # "miR-122" has 2 segments and would otherwise lose its family id.
    # "miR-16-1" has 3 segments -> paralog id "1", family "miR-16".
    # "let-7a-2" has 3 segments -> paralog id "2", family "let-7a".
    if len(segs) >= 3 and segs[-1].isdigit():
        paralog_id = segs[-1]
        family = "-".join(segs[:-1])
    else:
        paralog_id = ""
        family = stem
    return family, paralog_id, arm


def parse_mirbase(fa_path: Path) -> pd.DataFrame:
    """Parse miRBase mature.fa. Returns DataFrame with columns:
    species, mature_name, family, paralog_id, arm, mimat_id, seq."""
    rows = []
    cur = None
    seq = []
    # File is text despite the .gz extension (the previous download was raw text)
    with open(fa_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    rows.append((*cur, "".join(seq)))
                # Header: >hsa-miR-122-5p MIMAT0000421 Homo sapiens miR-122-5p
                parts = line[1:].split()
                full = parts[0]
                mimat = parts[1] if len(parts) > 1 else ""
                m = re.match(r"^([a-z]{3,4})-(.+)$", full)
                if not m:
                    cur = None
                    seq = []
                    continue
                species = m.group(1)
                rest = m.group(2)  # e.g., miR-122-5p, let-7a-2-3p
                family, paralog_id, arm = split_family_paralog_arm(rest)
                cur = (species, full, family, paralog_id, arm, mimat)
                seq = []
            else:
                seq.append(line.strip())
        if cur is not None:
            rows.append((*cur, "".join(seq)))
    df = pd.DataFrame(rows, columns=["species", "mature_name", "family",
                                     "paralog_id", "arm", "mimat_id", "seq"])
    return df


def build_orthologs(mb: pd.DataFrame) -> pd.DataFrame:
    """Pair hsa <-> mmu within (family, paralog_id, arm). Returns long-form pairs."""
    hsa = mb[mb["species"] == "hsa"].rename(columns={
        "mature_name": "human_mature", "mimat_id": "human_mimat",
        "seq": "human_seq"
    })[["family", "paralog_id", "arm",
        "human_mature", "human_mimat", "human_seq"]]
    mmu = mb[mb["species"] == "mmu"].rename(columns={
        "mature_name": "mouse_mature", "mimat_id": "mouse_mimat",
        "seq": "mouse_seq"
    })[["family", "paralog_id", "arm",
        "mouse_mature", "mouse_mimat", "mouse_seq"]]
    pairs = hsa.merge(mmu, on=["family", "paralog_id", "arm"], how="inner")
    return pairs


def gencode_mirna_symbols(meta_path: Path, biotype_col: str,
                          symbol_col: str, id_col: str) -> pd.DataFrame:
    """Return mapping of miRNA gene symbols (case-normalized) -> ENSG/ENSMUSG."""
    df = pd.read_csv(meta_path, sep=None, engine="python")
    df = df[df[biotype_col].astype(str) == "miRNA"].copy()
    df["symbol_norm"] = df[symbol_col].astype(str).str.upper()
    return df[["symbol_norm", symbol_col, id_col, biotype_col]].rename(
        columns={symbol_col: "_sym_raw", id_col: "_gene_id",
                 biotype_col: "_biotype"}
    )


def _stem_for_symbol(family: str) -> str:
    """Convert miRBase family stem to the alphanumeric core used by GENCODE
    symbols.

    Examples:
      "miR-122"   -> "122"            (HGNC: MIR122,    MGI: Mir122)
      "let-7a"    -> "let7a"          (HGNC: MIRLET7A,  MGI: Mirlet7a)
      "miR-30c"   -> "30c"            (HGNC: MIR30C,    MGI: Mir30c)
      "miR-100"   -> "100"
      "miR-let-7a"-> "let7a"          (safety: handle "let-7" prefixed with miR-)

    Returns the core lowercase string ready for case-flip into HGNC/MGI form.
    """
    f = family.lower().strip()
    # Strip leading "mir-" (the genus prefix)
    if f.startswith("mir-"):
        f = f[4:]
    # If what's left starts with "let-", drop the hyphen so "let-7a" -> "let7a"
    if f.startswith("let-"):
        f = "let" + f[4:]
    # Final hyphen scrub (e.g., "let-7a" could remain "let-7a" if no leading mir-)
    f = f.replace("-", "")
    return f


def normalize_mirna_to_gencode_symbol(family: str,
                                      paralog_id: str = "") -> list[str]:
    """miRBase family + paralog id -> list of candidate GENCODE symbols to try.

    Conventions (verified against GENCODE v49 / vM38):
      HGNC: paralog id is CONCATENATED, no hyphen.
            MIR122, MIRLET7A, MIR30A,
            MIRLET7A1 / MIRLET7A2 / MIRLET7A3,
            MIR30C1 / MIR30C2,
            MIR19B1 (some loci use the suffix, some don't -- emit both forms).
      MGI:  paralog id is HYPHENATED.
            Mir122, Mirlet7a, Mir30a,
            Mirlet7a-1 / Mirlet7a-2,
            Mir30c-1 / Mir30c-2.

    Strategy:
      - Always emit single-form (no paralog suffix) since paralogs share the
        same mature sequence at the family-arm level.
      - If paralog_id is present, ALSO emit cluster forms in both HGNC
        (concatenated) and MGI (hyphenated) styles, plus a couple of
        defensive variants (HGNC with hyphen, MGI concatenated) to catch
        any annotation drift.
    """
    stem = _stem_for_symbol(family)  # e.g., "122", "let7a", "30c"
    is_let = stem.startswith("let")
    hgnc_prefix = "MIRLET" if is_let else "MIR"
    mgi_prefix = "Mirlet" if is_let else "Mir"
    # The post-prefix core (e.g., "7a", "122", "30c")
    core = stem[3:] if is_let else stem
    hgnc_core = core.upper()
    mgi_core = core  # already lowercase
    candidates: set[str] = set()
    # Single-form (no paralog suffix) -- always emitted
    candidates.add(f"{hgnc_prefix}{hgnc_core}")
    candidates.add(f"{mgi_prefix}{mgi_core}")
    # Cluster-form variants
    if paralog_id:
        # HGNC canonical: concatenated, no hyphen (MIRLET7A2, MIR30C2)
        candidates.add(f"{hgnc_prefix}{hgnc_core}{paralog_id}")
        # MGI canonical: hyphenated (Mirlet7a-2, Mir30c-2)
        candidates.add(f"{mgi_prefix}{mgi_core}-{paralog_id}")
        # Defensive: HGNC with hyphen (some HGNC entries may use it)
        candidates.add(f"{hgnc_prefix}{hgnc_core}-{paralog_id}")
        # Defensive: MGI without hyphen (annotation drift)
        candidates.add(f"{mgi_prefix}{mgi_core}{paralog_id}")
    return sorted(candidates)


def _is_hgnc_form(sym: str) -> bool:
    """Heuristic: HGNC symbols are all-uppercase (MIR122, MIRLET7A2)."""
    return sym.isupper()


def _is_mgi_form(sym: str) -> bool:
    """Heuristic: MGI mouse symbols start uppercase then lowercase
    (Mir122, Mirlet7a, Mirlet7a-2)."""
    return (len(sym) >= 2 and sym[0].isupper()
            and sym[1:].lower() == sym[1:])


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mirbase", type=Path, default=DEFAULT_FA)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = p.parse_args()

    print(f"[mirbase] parsing {args.mirbase}", file=sys.stderr)
    mb = parse_mirbase(args.mirbase)
    print(f"[mirbase]  total miRNAs: {len(mb)}", file=sys.stderr)
    print(f"[mirbase]  hsa: {(mb['species']=='hsa').sum()}, "
          f"mmu: {(mb['species']=='mmu').sum()}", file=sys.stderr)

    pairs = build_orthologs(mb)
    print(f"[mirbase] hsa<->mmu pairs (family+paralog+arm matched): "
          f"{len(pairs)}", file=sys.stderr)

    # Map to GENCODE symbols (Ensembl IDs) where possible
    hsa_g = gencode_mirna_symbols(HUMAN_GTF_META, "gene_biotype",
                                  "gene_name", "ensembl_base")
    mmu_g = gencode_mirna_symbols(MOUSE_GTF_META, "mouse_biotype",
                                  "mouse_symbol_gtf", "mouse_ensembl_base")

    # Build candidate symbols per (family, paralog_id) and explode for matching.
    pairs["candidates"] = pairs.apply(
        lambda r: normalize_mirna_to_gencode_symbol(r["family"], r["paralog_id"]),
        axis=1,
    )
    # Explode into one row per candidate then match case-insensitively.
    long = pairs[["family", "paralog_id", "arm", "candidates"]].explode("candidates")
    long["symbol_norm"] = long["candidates"].astype(str).str.upper()

    # Human side: match symbol uppercase. Restrict candidates to HGNC-form
    # ones (all-uppercase) BEFORE the merge so we don't accidentally match a
    # mouse symbol that happens to live in a human-only metadata join (the
    # symbol_norm uppercase folds case anyway, so this is a sanity guard).
    hsa_g["symbol_norm"] = hsa_g["_sym_raw"].astype(str).str.upper()
    h_long = long[long["candidates"].apply(_is_hgnc_form)].copy()
    h_match = h_long.merge(
        hsa_g[["symbol_norm", "_sym_raw", "_gene_id"]].rename(
            columns={"_sym_raw": "human_symbol", "_gene_id": "human_ensembl"}),
        on="symbol_norm", how="inner"
    )[["family", "paralog_id", "arm", "human_symbol", "human_ensembl"]]
    # Mouse side: match symbol uppercase against MGI-form candidates.
    mmu_g["symbol_norm"] = mmu_g["_sym_raw"].astype(str).str.upper()
    m_long = long[long["candidates"].apply(_is_mgi_form)].copy()
    m_match = m_long.merge(
        mmu_g[["symbol_norm", "_sym_raw", "_gene_id"]].rename(
            columns={"_sym_raw": "mouse_symbol", "_gene_id": "mouse_ensembl"}),
        on="symbol_norm", how="inner"
    )[["family", "paralog_id", "arm", "mouse_symbol", "mouse_ensembl"]]

    # When we have a paralog id, prefer the cluster-form Ensembl match
    # (Mirlet7a-2) over the single-form (Mirlet7a) within that paralog row.
    # We do this by ranking cluster matches higher and keeping the first per
    # (family, paralog_id, arm) AFTER outer-joining the paralog rows. But to
    # preserve multiple Ensembl IDs across paralogs (Bug C fix), we DO NOT
    # dedup at the (family, arm) level any more -- we dedup at the per-row
    # input level only (each input (family, paralog_id, arm) row gets one
    # best human and one best mouse Ensembl match).
    def _rank_match_for_paralog(row, ensembl_col, symbol_col):
        # Higher rank = better. Cluster-form symbols (containing the paralog
        # id) outrank single-form when paralog_id is set.
        if not row["paralog_id"]:
            return 0
        sym = str(row[symbol_col]) if pd.notna(row[symbol_col]) else ""
        return 1 if row["paralog_id"] in sym else 0

    if not h_match.empty:
        h_match["_rank"] = h_match.apply(
            lambda r: _rank_match_for_paralog(r, "human_ensembl", "human_symbol"),
            axis=1,
        )
        h_match = (h_match.sort_values("_rank", ascending=False)
                          .drop_duplicates(["family", "paralog_id", "arm"],
                                           keep="first")
                          .drop(columns=["_rank"]))
    if not m_match.empty:
        m_match["_rank"] = m_match.apply(
            lambda r: _rank_match_for_paralog(r, "mouse_ensembl", "mouse_symbol"),
            axis=1,
        )
        m_match = (m_match.sort_values("_rank", ascending=False)
                          .drop_duplicates(["family", "paralog_id", "arm"],
                                           keep="first")
                          .drop(columns=["_rank"]))

    pairs = pairs.merge(h_match, on=["family", "paralog_id", "arm"], how="left")
    pairs = pairs.merge(m_match, on=["family", "paralog_id", "arm"], how="left")

    # NOTE on single-form rows whose family has paralogs (e.g., let-7a, where
    # HGNC has MIRLET7A1/2/3 but no bare MIRLET7A): we leave the Ensembl IDs
    # NaN on those rows. integrate_layers.py drops NaN-Ensembl rows, so the
    # single-form symbol (MIRLET7A) does not survive in the master table on
    # its own. ortholog_lookup.py handles single-form queries via a paralog-
    # prefix fallback (e.g., MIRLET7A -> match any MIRLET7A\d+ row).

    # Fallback symbols (use cluster-form candidate where available, else single).
    def _human_guess(row):
        cs = row["candidates"]
        # Prefer cluster-form HGNC if paralog_id present
        if row["paralog_id"]:
            for c in cs:
                if _is_hgnc_form(c) and row["paralog_id"] in c:
                    return c
        for c in cs:
            if _is_hgnc_form(c):
                return c
        return cs[0] if cs else ""

    def _mouse_guess(row):
        cs = row["candidates"]
        if row["paralog_id"]:
            for c in cs:
                if _is_mgi_form(c) and row["paralog_id"] in c:
                    return c
        for c in cs:
            if _is_mgi_form(c):
                return c
        return cs[0] if cs else ""

    pairs["human_gencode_symbol_guess"] = pairs.apply(_human_guess, axis=1)
    pairs["mouse_gencode_symbol_guess"] = pairs.apply(_mouse_guess, axis=1)

    print(f"[mirbase] pairs with both GENCODE Ensembl IDs: "
          f"{pairs[['human_ensembl','mouse_ensembl']].dropna().shape[0]}",
          file=sys.stderr)
    print(f"[mirbase] pairs with human-only Ensembl: "
          f"{(pairs['human_ensembl'].notna() & pairs['mouse_ensembl'].isna()).sum()}",
          file=sys.stderr)
    print(f"[mirbase] pairs with mouse-only Ensembl: "
          f"{(pairs['human_ensembl'].isna() & pairs['mouse_ensembl'].notna()).sum()}",
          file=sys.stderr)

    # Common-schema output -- one row per (family, paralog_id, arm) at this
    # point; we collapse arms below.
    out = pd.DataFrame({
        "mouse_ensembl": pairs["mouse_ensembl"],
        "mouse_symbol": pairs["mouse_symbol"].fillna(
            pairs["mouse_gencode_symbol_guess"]),
        "mouse_biotype": "miRNA",
        "human_ensembl": pairs["human_ensembl"],
        "human_symbol": pairs["human_symbol"].fillna(
            pairs["human_gencode_symbol_guess"]),
        "human_biotype": "miRNA",
        "tier_H_biomart": 0,
        "tier_H_mirbase": 1,
        "tier_H_mirgenedb": 0,
        "tier_M_phasej_synteny": 0,
        "tier_M_lncbook": 0,
        "tier_M_blast": 0,
        "tier_M_pseudogene_parent": 0,
        "tier_L_liftover": 0,
        "confidence_tier": "H",
        "evidence_count": 1,
        "provenance_sources": "mirbase_22_family_id",
        "mirbase_family": pairs["family"],
        "mirbase_paralog_id": pairs["paralog_id"],
        "mirbase_arm": pairs["arm"],
        "mirbase_human_mimat": pairs["human_mimat"],
        "mirbase_mouse_mimat": pairs["mouse_mimat"],
        "mirbase_human_mature_name": pairs["human_mature"],
        "mirbase_mouse_mature_name": pairs["mouse_mature"],
        "notes": "miRBase 22 family-ID matched",
    })

    # Collapse arms at the (mouse_ensembl, human_ensembl) gene-pair level.
    # Per requirement #4: for the integration with master table, the
    # (mouse_ensembl, human_ensembl) pair is the key; 5p and 3p of the same
    # gene pair collapse with mirbase_arm concatenated as "3p,5p".
    # Rows without BOTH ensembl IDs cannot be safely collapsed at the gene
    # level, so they are emitted as-is (one per arm). This matches the
    # downstream load_mirbase() behavior in integrate_layers.py which only
    # merges rows with both Ensembl IDs.
    key = ["mouse_ensembl", "human_ensembl"]
    both = out.dropna(subset=key, how="any").copy()
    missing = out[out["mouse_ensembl"].isna() | out["human_ensembl"].isna()].copy()

    def _csv_concat(series: pd.Series) -> str:
        vals = sorted({str(v) for v in series.dropna()
                       if str(v) != "" and str(v) != "nan"})
        return ",".join(vals)

    def _first_nonnull(series: pd.Series):
        for v in series:
            if pd.notna(v) and str(v) != "" and str(v) != "nan":
                return v
        return pd.NA

    concat_cols = ["mirbase_arm", "mirbase_paralog_id",
                   "mirbase_human_mimat", "mirbase_mouse_mimat",
                   "mirbase_human_mature_name", "mirbase_mouse_mature_name"]
    first_cols = ["mouse_symbol", "mouse_biotype", "human_symbol",
                  "human_biotype", "tier_H_biomart", "tier_H_mirbase",
                  "tier_H_mirgenedb", "tier_M_phasej_synteny",
                  "tier_M_lncbook", "tier_M_blast",
                  "tier_M_pseudogene_parent", "tier_L_liftover",
                  "confidence_tier", "evidence_count",
                  "provenance_sources", "mirbase_family", "notes"]
    agg_dict = {c: _csv_concat for c in concat_cols if c in both.columns}
    agg_dict.update({c: _first_nonnull for c in first_cols if c in both.columns})
    if not both.empty:
        collapsed = both.groupby(key, as_index=False).agg(agg_dict)
    else:
        collapsed = both

    final = pd.concat([collapsed, missing], ignore_index=True, sort=False)
    # Reorder columns to match the original schema (plus mirbase_paralog_id).
    col_order = [
        "mouse_ensembl", "mouse_symbol", "mouse_biotype",
        "human_ensembl", "human_symbol", "human_biotype",
        "tier_H_biomart", "tier_H_mirbase", "tier_H_mirgenedb",
        "tier_M_phasej_synteny", "tier_M_lncbook", "tier_M_blast",
        "tier_M_pseudogene_parent", "tier_L_liftover",
        "confidence_tier", "evidence_count", "provenance_sources",
        "mirbase_family", "mirbase_paralog_id", "mirbase_arm",
        "mirbase_human_mimat", "mirbase_mouse_mimat",
        "mirbase_human_mature_name", "mirbase_mouse_mature_name",
        "notes",
    ]
    col_order = [c for c in col_order if c in final.columns]
    final = final[col_order]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(args.out, sep="\t", index=False)
    print(f"[mirbase] wrote {len(final)} miRNA pairs (after arm collapse) "
          f"to {args.out}", file=sys.stderr)
    print(f"[mirbase]  with both Ensembl IDs: "
          f"{final[['mouse_ensembl','human_ensembl']].dropna().shape[0]}",
          file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
